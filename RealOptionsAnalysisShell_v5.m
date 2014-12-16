%-------------------------------------------------------------------------%
%-----------------------Real Options Analysis Shell-v4--------------------%
%-------------------------------------------------------------------------%
%
%This code performs a Real Options valuation of an investment scenario
%following the proceedures outlined in Purvis Et Al 1995 with additional
%code to simulate an industry-wide adoption model with learning spillovers.
%
%Removed commented and defunct code from v3
%Removed large dupicate code under SIM section: if SIM = 1, else
%Originally named ReturnsToInvestment_post_PurvisEtAl_v5
%
%-------------------------------------------------------------------------%
%---Written by Gregory M. Astill, Nov. 14, 2014, for Matlab---------------%
%-------------------------------------------------------------------------%



%-----------------------Preamble------------------------------------------%

clear
clc
format compact

graphs = 0;                     %Turns on/off graph (mc_sigma)



%-----------------------Model Varaibles-----------------------------------%

N = 100;                        %Number of firms
rho = .06                       %Discount Rate
MC= 200;                        %Number of Monte Carlo simulations in mc_sigma
N1= 15;                         %Investment lifetime
J   = 4;                        %Number of firms in the industry
subsidy = 0;                    %On/Off for government grants
weig = .5;                      %Arbitrary weight
socben = subsidy*(weig*3.4566+(1-weig)*6.9132)
                                %Arbitrary calculation of social benefit
ejbar = 0;                      %Proportion of firms that have adopted
SIM = 6;                        %Number of years of firm adoption simulation
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
NnextState = Dataout(:,4);      %One period number of AD owners in state next door draw
EnvR = Dataout(:,5);            %One period environmental evaluations of AD benefits draw
sEnvR = var(EnvR);              %Spread of AD evaluation among all periods         
u = round(100*unifrnd(0.01,N/100,J,1));
firmEnvR = EnvR(u,1)            %Firm constant environmental valuation

expR = sum(R)/N                 %Expected revenue
Kout = Kout - socben;           %Capital cost subsidized by the amount of social benefit
K = sum(Kout)/N;                %Average K
M = rho*K                       %Marshellian trigger for infinite lifetime
Mstar= triggersolve(rho,N1,K)   %Marshellian trigger for N1 year lifetime
              


%-----------------------Analysis------------------------------------------%
%-----------------------Analysis: Loop over learning from others and no learning----%

for al=0:1               %%%%%%%% TO DO, insert another for loop for subsidies so all EnvR are the same 
alpha=al;                       %Learning from others if alpha = 1
firmj = zeros(J,5,SIM+1);
trialmat = zeros(J*SIM,14);     %Output data

%-----------------------Analysis: Perform SIM periods of adoption model%

    for sim=1:SIM

        d = ones(J*SIM,1);
        trialmat(:,1) = (trial(1,1)+1);          %Number of time this code was run
        uncertain = var(Kout);                  %This is defined on line 29 in mc_sigma_v7.m
        trialmat(:,[3:6]) = d*[K, expR, rho, uncertain];
        trialmat(:,14) = sEnvR;               %Spread of AD evaluation among all periods
        trialmat(:,15) = socben;                %Social benefit or equivalent subsidy
        trialmat(:,16) = alpha;                 %Learning from others
        trialmat([1+(sim-1)*J:sim*J],2) = sim;        %Record sim period
%          trialmat(:,2)
%          pause
        for j=1:J
            trialmat((sim-1)*J+j,7) = j;          %Record firm
            trialmat((sim-1)*J+j,13) = firmEnvR(j,1); %Record firm's environmental valuation
                
            if firmj(j,1,sim)==0                %Firms haven't adopted
                REnvR = R + firmEnvR(j,1)
                s = mc_sigma_v10(REnvR, rho, ejbar, alpha, N1, K, Kout, N, MC, graphs);            
                                                %Indexing is off between firmj and trialmat
                                                %due to the need to reference 
                firmj(j,2,sim+1) = s;           %last years' adoption status
                trialmat((sim-1)*J+j,8) = s;    %Variance term sigma of Brownian motion

                beta = (1/2)*(1+sqrt(1+(8*rho)/s));
                firmj(j,3,sim+1) = beta;          %Beta term of modified investment trigger
                trialmat((sim-1)*J+j,9) = beta;

                rhop = ((beta/(beta - 1))*rho);
                firmj(j,4,sim+1) = rhop;          %Modified hurdle discount rate
                trialmat((sim-1)*J+j,10) = rhop;

                H = rhop*K;                     %Average modified investment trigger
                Hstar = triggersolve(rhop,N1,K) %Modified investment trigger
                firmj(j,5,sim+1) = Hstar;
                trialmat((sim-1)*J+j,11) = Hstar;
                pause
                if Hstar < expR + firmEnvR(j,1)     %Firm j's adoption decision

                    firmj(j,1,sim+1) = 1;         
                    trialmat((sim-1)*J+j,12) = 1;
                    disp('ej = 1')
                end

            else                                %Firms have already adopted

                trialmat((sim-1)*J+j,8) = -99;
                trialmat((sim-1)*J+j,9) = -99;
                trialmat((sim-1)*J+j,10) =-99;
                trialmat((sim-1)*J+j,11) =-99;
                firmj(j,1,sim+1) = 1;
                trialmat((sim-1)*J+j,12) = 1;
            end           

        end

        ejbar = mean(firmj(:,1,sim+1))          %Calculate new adoption ratio
        disp(' ')
        disp([sim, ' simulations'])
    
    end

AA = sprintf('A%d', alpha*(SIM*J) + trial(1,2) + 1)
OO = sprintf('P%d', alpha*(SIM*J) + SIM*J + trial(1,2))
trialrange = sprintf('%s:%s', AA, OO)
%xlswrite('C:\Users\Greg\Google Drive\WSU_F12\ExAnteTechAdoptionUnderUncertainty\Fall2014RevisionUsingCowleyData\trials1.xlsx',...
%    trialmat,trialrange)
xlswrite(resultsfile,trialmat,trialrange)
%SIM = SIMM;
out = firmj;
out

%print(fig,filename,'-dpng','-r300')

end







