%-------------------------------------------------------------------------%
%-----------------------ROASimulator_v1-----------------------------------%
%-------------------------------------------------------------------------%
%
%This code calls the RealOptionsAnalysisShell_v10 function to perform ROA
%under various parameter settings designated herein.
%
%-------------------------------------------------------------------------%
%---Written by Gregory M. Astill, Nov. 14, 2014, for Matlab---------------%
%-------------------------------------------------------------------------%

graphonoff = 0;
N = 100;                        %Number of distribution draws from data
rho = .04;                       %Discount Rate
MC= 2500;                        %Number of Monte Carlo simulations in mc_sigma
N1= 15;                         %Investment lifetime
J   = 50;                        %Number of firms in the industry
subsidy = 0;                    %On/Off for government grants
eps = 3;                        %Variance adjustor from uncertainty re: 
                                %revenue distribution. Very small eps,
                                %increases variance. Very large, decreases
SIM = 15;                        %Number of years of firm adoption simulation
aK = 2 ;                        %Adjusts K   -- 1 means no adjustment
aE = 00;                        %Adjusts EnvR-- 1 means no adjustment

%output=zeros(9,5,1,1);

for irho=1:9
rho = irho/100;
    for iN1=3:6
    N1 = iN1*5;
        for ieps=1:10
        eps = ieps/2;
            for iaK=0:10
            aK = iaK/4;
                for iaE=0:10
                aE = iaE/10;
    %trialmat = RealOptionsAnalysisShell_v9(rho)
    % trialmat = RealOptionsAnalysis_v10(rho);
    trialmat = RealOptionsAnalysis_v10(graphonoff, N, rho, MC, N1, ...
                            J, subsidy, eps, SIM, aK, aE);
    %[r,c] = size(trialmat);
    %output(irho,iN1,:,:) = trialmat;
                end
            end
        end
    end
end
