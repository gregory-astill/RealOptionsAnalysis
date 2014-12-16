%-------------------------------------------------------------------------%
%-----------------------Real Options Analysis Shell-v4--------------------%
%-------------------------------------------------------------------------%
%
%Removed commented and defunct code from v3
%Originally named ReturnsToInvestment_post_PurvisEtAl_v5
%
%This code performs a Real Options valuation of an investment scenario
%following the proceedures outlined in Purvis Et Al 1995 with additional
%code to simulate an industry-wide adoption model with learning spillovers.
%
%-------------------------------------------------------------------------%
%---Written by Gregory M. Astill, Nov. 14, 2014, for Matlab---------------%
%-------------------------------------------------------------------------%



%-----------------------Preamble------------------------------------------%

clear
clc
format compact

graphs = 0;



%-----------------------Model Varaibles-----------------------------------%

N = 100;                        %Number of firms
rho = .03                       %Discount Rate
MC= 200;                        %Number of Monte Carlo simulations in mc_sigma
N1= 15;                         %Investment lifetime
J   = 2;                        %Number of firms in the industry
subsidy = 0;                    %On/Off for government grants
weig = .5;                      %Arbitrary weight
socben = subsidy*(weig*3.4566+(1-weig)*6.9132)
                                %Arbitrary calculation of social benefit
SIM = 1;                        %Number of Monte Carlo Simulations
SIMM = 3;                       %After the first Monte Carlo simulation by SIM
                                %increase thereafter to SIMM
Hstar = 0;                      %Set initial value of modified trigger


%-----------------------Data----------------------------------------------%

resultsfile = 'C:\Users\Greg\Google Drive\WSU_F12\ExAnteTechAdoptionUnderUncertainty\Fall2014RevisionUsingCowleyData\trials1.xlsx';
trial = xlsread(resultsfile, 'B1:C1') 
                                %Appends results from a single run of this code into an excel workbook.

datafile = 'C:\Users\Greg\Google Drive\ROA&LearningSpillovers\CowleyData\CowleyDataOut.csv';
datafile2 = 'C:\Users\Greg\Google Drive\ROA&LearningSpillovers\CowleyData\CowleyDataOut_HeterogeneityOfEnvironmentalValuation.csv';


Dataout = DataPrep(datafile, datafile2, N, subsidy);
                                %Reads in distribution data from two tables and make N draws randomly.

Kout = Dataout(:,1);            %One period capital draw
R = Dataout(:,2);               %One period revenue draw
NinState = Dataout(:,3);        %One period number of AD owners in state draw
NnextState = Dataout(:,4);      %One period number of AD owners in state next door 
EnvR = Dataout(:,5);            %One period environmental evaluations of AD benefits
sEnvR = var(EnvR);              %Spread of AD evaluation among all periods

expR = sum(R)/N                 %Expected revenue
Kout = Kout - socben;           %Capital cost subsidized by the amount of social benefit
K = sum(Kout)/N;                %Average K
M = rho*K                       %Marshellian trigger for infinite lifetime
Mstar= triggersolve(rho,N1,K)   %Marshellian trigger for N1 year lifetime
              


%-----------------------Analysis------------------------------------------%
%-----------------------Analysis: Loop over subsidies and no subsidies----%

for al=0:1          %%%%%%%% TO DO, insert another for loop for subsidies so all EnvR are the same 
alpha=al;        %Learning from others if alpha = 1
firmj = zeros(J,5,SIM);


%-----------------------Analysis: Perform SIM Monte Carlos of adoption model%
for sim=1:SIM
    if sim == 1
    ejbar = 0;                  % Proportion of firms that have adopted
        trialmat = zeros(J*SIM,14);             %for output file
        trialmat(:,1) = trial(1,1)+1;
        for i=1:SIM
            trialmat([1+(i-1)*J:i*J],2) = i;
            for j=1:J
                trialmat((i-1)*J+j,7) = j;
                trialmat((i-1)*J+j,13) = EnvR(j,1); 
            end
        end
        uncertain = var(Kout); %This is defined on line 29 in mc_sigma_v7.m
        d = ones(J*SIM,1);
        trialmat(:,[3:6]) = d*[K, expR, rho, uncertain];
        trialmat(:,14) = d*sEnvR;
        trialmat(:,15) = socben;

    
    for j=1:J
        
        
        REnvR = R + EnvR;
        s = mc_sigma_v10(REnvR, rho, ejbar, alpha, N1, K, Kout, N, MC);            


        firmj(j,2,sim) = s;
        trialmat((sim-1)*J+j,8) = s;
        %s = 0.043
        beta = (1/2)*(1+sqrt(1+(8*rho)/s));
        firmj(j,3,sim) = beta/(beta-1);
        trialmat((sim-1)*J+j,9) = beta/(beta-1);
        rhop = ((beta/(beta - 1))*rho);
        firmj(j,4,sim) = rhop;
        trialmat((sim-1)*J+j,10) = rhop;
        H = rhop*K;
        Hstar = triggersolve(rhop,N1,K)
        firmj(j,5,sim) = Hstar;
        trialmat((sim-1)*J+j,11) = Hstar;

        if Hstar < expR + EnvR(j,1)

            firmj(j,1,sim) = 1;
            trialmat((sim-1)*J+j,12) = 1;
            disp('ej = 1')
        end

        
        

    end
    
    else
    ejbar = mean(firmj(:,1,sim-1))
        
    for j=1:J
        jth = j
        if firmj(j,1,sim-1)==0
            
            REnvR = R + EnvR(j,1);
            s = mc_sigma_v10(REnvR, rho, ejbar, alpha, N1, K, Kout, N, MC);
   
            
            firmj(j,2,sim) = s;
            trialmat((sim-1)*J+j,8) = s;
            %s = 0.043
            beta = (1/2)*(1+sqrt(1+(8*rho)/s));
            firmj(j,3,sim) = beta/(beta-1);
            trialmat((sim-1)*J+j,9) = beta/(beta-1);
            rhop = ((beta/(beta - 1))*rho);
            firmj(j,4,sim) = rhop;
            trialmat((sim-1)*J+j,10) = rhop;
            H = rhop*K;
            Hstar = triggersolve(rhop,N1,K)
            firmj(j,5,sim) = Hstar;
            trialmat((sim-1)*J+j,11) = Hstar;
            
           
            if Hstar < expR + EnvR(j,1)
                
                firmj(j,1,sim) = 1;
                trialmat((sim-1)*J+j,12) = 1;
                disp('ej = 1')
            end
            
            
        else
            firmj(j,1,sim) = 1;
            trialmat((sim-1)*J+j,12) = 1;
            
        end
    end
    
    end
    
    disp(' ')
    disp([sim, ' simulations'])
    
end

AA = sprintf('A%d', (1-alpha)*(trial(1,2)+1) + alpha*(trial(1,2)+J+1))
OO = sprintf('O%d', (1-alpha)*(trial(1,2)+J) + alpha*(trial(1,2)+J+(J*SIM)))
trialrange = sprintf('%s:%s', AA, OO)
xlswrite('C:\Users\Greg\Google Drive\WSU_F12\ExAnteTechAdoptionUnderUncertainty\Fall2014RevisionUsingCowleyData\trials1.xlsx',...
    trialmat,trialrange)
SIM = SIMM;
out = firmj;
out

end







