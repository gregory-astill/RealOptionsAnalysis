    %v5 begun 11/11/14
clear
clc
format compact

%ADbata = csvread('revenue_only.csv', 0, 0, [0, 0, 4, 3]);
trial = xlsread('C:\Users\Greg\Google Drive\WSU_F12\ExAnteTechAdoptionUnderUncertainty\Fall2014RevisionUsingCowleyData\trials1.xlsx', 'B1:C1')

rho = .03
SIMM = 3;

w1 = 1;
eps= w1*normrnd(0,1);
N1= 15;     %number of years the investment will give returns
J   = 20;
% muEnvR = 1;
% sEnvR = .5; %added below
% EnvR = normrnd(muEnvR,sEnvR,1,J); %arbitrarily chosen
% EnvR = sort(EnvR,'descend');
% EnvR = 0;
%EnvR = zeros(1,J);
%min = min(EnvR)
subsidy = 1;
weig = .5;
socben = subsidy*(weig*3.4566+(1-weig)*6.9132)


subsidy = 00;

%K1 = ADbata(:,3) 
%K  = 996.2

N = 100;
datafile = 'C:\Users\Greg\Google Drive\ROA&LearningSpillovers\CowleyData\CowleyDataOut.csv';
datafile2 = 'C:\Users\Greg\Google Drive\ROA&LearningSpillovers\CowleyData\CowleyDataOut_HeterogeneityOfEnvironmentalValuation.csv';
Dataout = DataPrep(datafile, datafile2, N);
%load Ksim.mat
Kout = Dataout(:,1);
R = Dataout(:,2);
NinState = Dataout(:,3);
NnextState = Dataout(:,4);
EnvR = Dataout(:,5);
sEnvR = var(EnvR);

eta2 = numel(Kout);
Kout = Kout - socben;
K = sum(Kout)/eta2
M = rho*K                   %Marshellian trigger = discount rate * capital investment
Mstar= triggersolve(rho,N1,K)
Hstar = 0;
%ER1 = ADbata(:,4)
%ER  = sum(ER1)/numel(ER1)   %Expected return
%R = ADbata(:,2);
%load Rsim.mat

%load Psim.mat
%N = 100;
a = 3;
% R = bsR(:,1);
%R = Rsim
%pause
[f,ri] = ksdensity(R);
plot(ri,f)
%pause

eta = numel(R)

expR = sum(R)/eta
SIM = 1;
for al=0:1          %%%%%%%% TO DO, insert another for loop for subsidies so all EnvR are the same 
alpha=al;        %Learning from others if alpha = 1
firmj = zeros(J,5,SIM);

for sim=1:SIM
    if sim == 1
    ejbar = 0;
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
        
        
        REnvR = R + EnvR(j,1);
        s = mc_sigma_v8(REnvR, rho, eta, ejbar, alpha, N1, eps, K, Kout, eta2);
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
            s = mc_sigma_v8(REnvR, rho, eta, ejbar, alpha, N1, eps, K, Kout, eta2);
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

%bo1mbeta = beta/(beta-1);
%out = [N1+1 H s rhop bo1mbeta];

AA = sprintf('A%d', (1-alpha)*(trial(1,2)+1) + alpha*(trial(1,2)+J+1))
OO = sprintf('O%d', (1-alpha)*(trial(1,2)+J) + alpha*(trial(1,2)+J+(J*SIM)))
trialrange = sprintf('%s:%s', AA, OO)
xlswrite('C:\Users\Greg\Google Drive\WSU_F12\ExAnteTechAdoptionUnderUncertainty\trials1.xlsx',...
    trialmat,trialrange)
SIM = SIMM;
out = firmj;
out

end







