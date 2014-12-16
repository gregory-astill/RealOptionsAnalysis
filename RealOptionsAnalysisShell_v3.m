    %v5 begun 11/11/14
    %Renamed from ReturnsToInvestment_post_PurvisEtAl_v5
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
J   = 2;   %numbre of firms in the industry
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
%[f,ri] = ksdensity(R);
%plot(ri,f)
%pause

eta = numel(R)

expR = sum(R)/eta
SIM = 1;
for al=0:1          %%%%%%%% TO DO, insert another for loop for subsidies so all EnvR are the same 
alpha=al;        %Learning from others if alpha = 1
firmj = zeros(J,5,SIM);

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
        %s = mc_sigma_v9(REnvR, rho, eta, ejbar, alpha, N1, eps, K, Kout, eta2);

%%----------------------------------------------%%       
%%Begin mc_sigma function-----------------------%%
%%----------------------------------------------%%
X = REnvR;
MC= 2000;
Np1=N1+1;


PVt = zeros(N1,1);
PVtp1= zeros(N1,1);
muv = zeros(MC,1);
expmuv= zeros(MC,1);
mcdeltaLnVjL = 0;
mcdeltaLnVj = zeros(MC,1);

sig2v = zeros(MC,1);
mcsig2vL = 0;
mcsig2v = zeros(MC,1);
expsig2v = zeros(MC,1);

DELTALnVj = zeros(MC,1);
jj=1;
errors = 0;
mc = 0;

    while mc < MC
        errori = 0;
        errorj = 0;
        limitij= 100;
        while jj==mc
        TX = zeros(Np1,1);
        tX = zeros(N1,1);
        tp1X=zeros(N1,1);
        u = round(100*unifrnd(0.01,eta/100,Np1,1));


        for n1=1:Np1   
        
            TX(n1,1) = X(u(n1,1));
            %blah = TX(n1,1)
            if n1 < Np1
                tX(n1,1) = TX(n1,1);
                PVt(n1,1) = tX(n1,1)/((1+rho)^(n1));
            end
            
            if n1 > 1
                tp1X(n1-1,1) = TX(n1,1);
                PVtp1(n1-1,1)=tp1X(n1-1,1)/((1+rho)^(n1-1));
            end
            
        end
        %pause
        
        

        if sum(PVt) < K
            errori = errori + 1;
            
        elseif sum(PVtp1) < K
            errorj = errorj + 1;
            
        else
        %comb1=[tX PVt tp1X PVtp1];
        omega0 = rho/(1-(1/((1+rho)^(N1-0))));      %t = 0 b/c reference is invest now
        omega1 = rho/(1-(1/((1+rho)^(N1-0-1))));    %t-1 
        Vt = log(omega0*sum(PVt)/rho)
        Vtp1=log(omega1*sum(PVtp1)/rho)           %updated in v7 to include learning from others
        deltaLnVj=Vt-Vtp1                         %   via investment cost uncertainty

        eps2 = normrnd(0,abs(deltaLnVj)/3);  %corresponds to Purvis et al. extension
                                             %the reduction factor that dimishes as more
                                             %firms adopt and others learn from them
        Vnoise = (1-ejbar*alpha)*eps2;             %    learning from others term that goes to zero as more firms adopt
        deltaLnVj = deltaLnVj + Vnoise
                            
%         if imag(Vt) ~= 0 | imag(Vtp1) ~= 0
%             warning('Imaginary numbers')
%             pause
%         end
    
       
        %{
        if mc < 3
            deltaLnVj
            eps
            ejbar
            yzxsa = eps*(1-ejbar)*alpha
        end
        %}

        DELTALnVj(mc,1) = deltaLnVj;
        %muv(mc,1) = deltaLnVj + eps*(1-ejbar)*alpha;    %Learning from others
        %mcdeltaLnVjL=mcdeltaLnVjL + deltaLnVj;
        %mcdeltaLnVj(mc,1)=mcdeltaLnVjL;
        %expmuv(mc,1)=(mcdeltaLnVj(mc,1)/mc);
        expmuv(mc,1)=sum(DELTALnVj)/mc;
        jj = jj + 1;
        end
        
        if errori > limitij | errorj > limitij
            display Error sum(PVt) < K or sum(PVtp1) < K more than 100 times.
            errors = 1;
            jj = jj + 1;
            mc = MC;
        end
       
%         TX
%         bl = sum(PVt)
%         blp1=sum(PVtp1)
%         Vt
%         Vtp1
%         DELTALnVj(mc,1)
%         expmuv(mc,1)     
%         pause
       
        end

    mc = mc + 1;   
    end

    for mc=1:MC
        
        sig2v(mc,1) = (DELTALnVj(mc,1) - expmuv(MC,1))^2;
%         mcsig2vL = mcsig2vL + sig2v(mc,1);
%         mcsig2v(mc,1) = mcsig2vL;
%         expsig2v(mc,1)=(mcsig2v(mc,1)/mc);
        expsig2v(mc,1)=sum(sig2v)/mc;
        
    end
   
expmu = expmuv(MC,1)
if errors == 0
s = expsig2v(MC,1)
else
s = 1
end

% count1=timeseries(expmuv(:,1),1:MC);
% count2=timeseries(expsig2v(:,1),1:MC);
% plot(count1, '-b')
% hold on
% plot(count2, '-m')
% title('Monte Carlo Simulation')
% xlabel('Monte Carlo Iterations')
% legend('mean','variance')


%%----------------------------------------------%%
%%--End mc_sigma function-----------------------%%
%%----------------------------------------------%%

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
            s = mc_sigma_v9(REnvR, rho, eta, ejbar, alpha, N1, eps, K, Kout, eta2);


            
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
xlswrite('C:\Users\Greg\Google Drive\WSU_F12\ExAnteTechAdoptionUnderUncertainty\Fall2014RevisionUsingCowleyData\trials1.xlsx',...
    trialmat,trialrange)
SIM = SIMM;
out = firmj;
out

end







