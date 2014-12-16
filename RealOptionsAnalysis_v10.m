function [ trialmat ] = RealOptionsAnalysis_v10(graphonoff, N, rho, MC, N1, ...
                        J, subsidy, eps, SIM, aK, aE)
                    
%-------------------------------------------------------------------------%
%-----------------------Real Options Analysis v10-------------------------%
%-------------------------------------------------------------------------%
%
%This code performs a Real Options valuation of an investment scenario
%following the proceedures outlined in Purvis Et Al 1995 with additional
%code to simulate an industry-wide adoption model with learning spillovers.
%
%Removed commented and defunct code from v3
%Removed large dupicate code under SIM section: if SIM = 1, else
%Originally named ReturnsToInvestment_post_PurvisEtAl_v5
%v6 allows adjustment to Cowley data as robustness checks
%v7 TODO runs learning from others on exact same data as no learning
%v8 removed firmj variable b/c it is redundant with trialmat
%-------------------------------------------------------------------------%
%---Written by Gregory M. Astill, Nov. 14, 2014, for Matlab---------------%
%-------------------------------------------------------------------------%

% graphonoff = 0;
% N = 100;                        %Number of distribution draws from data
% rho = .04                       %Discount Rate
% MC= 150;                        %Number of Monte Carlo simulations in mc_sigma
% N1= 15;                         %Investment lifetime
% J   = 3;                        %Number of firms in the industry
% subsidy = 0;                    %On/Off for government grants
weig = .5;                      %Arbitrary weight
socben = subsidy*(weig*3.4566+(1-weig)*6.9132);
                                %Arbitrary calculation of social benefit
ejbar = 0;                      %Proportion of firms that have adopted
% eps = 3;                        %Variance adjustor from uncertainty re: 
                                %revenue distribution. Very small eps,
                                %increases variance. Very large, decreases
% SIM = 2;                        %Number of years of firm adoption simulation
               
Hstar = 0;                      %Set initial value of modified trigger

% aK = 2 ;                        %Adjusts K   -- 1 means no adjustment
% aE = 00;                        %Adjusts EnvR-- 1 means no adjustment

%-----------------------Data----------------------------------------------%

%Manage results in excel
resultsfile = 'C:\Users\Greg\Google Drive\WSU_F12\ExAnteTechAdoptionUnderUncertainty\Fall2014RevisionUsingCowleyData\trials1.xlsx';
resultsopen = 0;                %Records whether trials1 is open in excel
try 
    Excel = actxGetRunningServer('Excel.Application');
    Workbooks = Excel.Workbooks;
        for ii = 1:Workbooks.Count
            if strcmp(resultsfile,Workbooks.Item(ii).FullName)
                resultsopen = 1;    %If trials1 is open, close it using xls_check_if_open
                                    %and reopen it at the bottom of this file
            end
        end
    xls_check_if_open(resultsfile,'close');
catch MExc
    Excel = actxGetRunningServer('Excel.Application');
    Workbooks = Excel.Workbooks;
end
    
trial = xlsread(resultsfile, 'B1:C1');
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
EnvR = aE .* EnvR;              %Adjusts EnvR 
sEnvR = var(EnvR);              %Spread of AD evaluation among all periods 
u = round(100*unifrnd(0.01,N/100,J,1));
firmEnvR = EnvR(u,1);            %Firm constant environmental valuation

expR = sum(R)/N;                 %Expected revenue
Kout = Kout - socben;           %Capital cost subsidized by the amount of social benefit
K = sum(Kout)/N;                %Average K
K = aK .* K;                    %Adjusted K
M = rho*K;                       %Marshellian trigger for infinite lifetime
Mstar= triggersolve(rho,N1,K);   %Marshellian trigger for N1 year lifetime
              


%-----------------------Analysis------------------------------------------%
%-----------------------Analysis: Loop over learning from others and no learning----%

for al=0:1               %%%%%%%% TO DO, insert another for loop for subsidies so all EnvR are the same 
alpha=al;                       %Learning from others if alpha = 1
ejbar = 0;                      %No one has adopted yet
trialmat = zeros(J*SIM,23);     %Output data

%-----------------------Analysis: Perform SIM periods of adoption model%

    for sim=1:SIM
        
        d = ones(J*SIM,1);
        trialmat(:,1) = (trial(1,1)+1);          %Number of time this code was run
        uncertain = var(Kout);                  %This is defined on line 29 in mc_sigma_v7.m
        trialmat(:,3:6) = d*[K, expR, rho, uncertain];
        trialmat(:,14) = sEnvR;               %Spread of AD evaluation among all periods
        trialmat(:,15) = socben;                %Social benefit or equivalent subsidy
        trialmat(:,16) = alpha;                 %Learning from others
        trialmat(:,17) = MC;                    %Number of Monte Carlos
        trialmat(:,18) = aE;
        trialmat(:,19) = aK;
        trialmat(:,20) = N1;
        trialmat(:,21) = N;
        trialmat(:,22) = subsidy;
        trialmat(:,23) = eps;
        trialmat(1+(sim-1)*J:sim*J,2) = sim;        %Record sim period
%          trialmat(:,2)
%          pause;
       
        
        for j=1:J
            trialmat((sim-1)*J+j,7) = j;          %Record firm
            trialmat((sim-1)*J+j,13) = firmEnvR(j,1); %Record firm's environmental valuation
                 
            if (sim-2)*J+j < 1 || trialmat((sim-2)*J+j,12) == 0 
                                                %Ignore first round because no one has adopted
                                                %Or after, choose Firms who haven't adopted                                
                trialmat((sim-1)*J+j,12);
                REnvR = R + firmEnvR(j,1);
                testR = normrnd(50,14,100,1);
                s = mc_sigma_v11(testR, rho, ejbar, alpha, N1, N, MC, graphonoff, eps);              
                %s = mc_sigma_v10(REnvR, rho, ejbar, alpha, N1, N, MC, graphonoff);            
                                                %Indexing is off between firmj and trialmat
                                                %due to the need to reference 
                                                %last years' adoption status
                trialmat((sim-1)*J+j,8) = s;    %Variance term sigma of Brownian motion

                beta = (1/2)*(1+sqrt(1+(8*rho)/s));
                trialmat((sim-1)*J+j,9) = beta; %Beta term of modified investment trigger
                
                rhop = ((beta/(beta - 1))*rho);
                trialmat((sim-1)*J+j,10) = rhop;%Modified hurdle discount rate

                H = rhop*K;                     %Average modified investment trigger
                Hstar = triggersolve(rhop,N1,K); %Modified investment trigger
                trialmat((sim-1)*J+j,11) = Hstar;
%                pause;
                if Hstar < expR + firmEnvR(j,1)     %Firm j's adoption decision        
                    trialmat((sim-1)*J+j,12) = 1;
%                     disp('ej = 1');
                end

            else                                %Firms have already adopted

                trialmat((sim-1)*J+j,8) = -99;
                trialmat((sim-1)*J+j,9) = -99;
                trialmat((sim-1)*J+j,10) =-99;
                trialmat((sim-1)*J+j,11) =-99;
                trialmat((sim-1)*J+j,12) = 1;
            end           

        end
        
%         sim
%         lowblah = (sim-1)*J+1
%         highblah= (sim-1)*J+J
        
                                               %Calculate new adoption ratio
        ejbar = mean(trialmat((sim-1)*J+1:(sim-1)*J+J,12));
%         disp(' ');
%         disp([sim, ' simulations']);
        
    
    end

AA = sprintf('A%d', alpha*(SIM*J) + trial(1,2) + 1);
OO = sprintf('W%d', alpha*(SIM*J) + SIM*J + trial(1,2));
trialrange = sprintf('%s:%s', AA, OO);
%xlswrite('C:\Users\Greg\Google Drive\WSU_F12\ExAnteTechAdoptionUnderUncertainty\Fall2014RevisionUsingCowleyData\trials1.xlsx',...
%    trialmat,trialrange)
xlswrite(resultsfile,trialmat,trialrange);

%out

    if resultsopen == 1                         %Reopen trials1 in Excel
        winopen(resultsfile);
    end

%print(fig,filename,'-dpng','-r300')

end

