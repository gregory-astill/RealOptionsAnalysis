function s = mc_sigma_v9(X, rho, ejbar, alpha, N1, K, Kout, N) 
    %v5 begun 5/24

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
j=1;
errors = 0;
mc = 0;

    while mc < MC
        errori = 0;
        errorj = 0;
        limitij= 100;
        while j==mc
        TX = zeros(Np1,1);
        tX = zeros(N1,1);
        tp1X=zeros(N1,1);
        u = round(100*unifrnd(0.01,N/100,Np1,1));
        %eps2 = normrnd(0,250);   %corresponds to Purvis et al. extension
        u2 = round(unifrnd(1,N));
        eps2 = Kout(u2)-K;

        for n1=1:Np1   
        
            TX(n1,1) = X(u(n1,1));
            
            if n1 < Np1
                tX(n1,1) = TX(n1,1);
                PVt(n1,1) = tX(n1,1)/((1+rho)^(n1));
            end
            
            if n1 > 1
                tp1X(n1-1,1) = TX(n1,1);
                PVtp1(n1-1,1)=tp1X(n1-1,1)/((1+rho)^(n1-1));
            end
            
        end
        
        

        if sum(PVt) < K
            errori = errori + 1;
            
        elseif sum(PVtp1) < K
            errorj = errorj + 1;
            
        else
        %comb1=[tX PVt tp1X PVtp1];
        omega0 = rho/(1-(1/((1+rho)^(N1-0))));      %t = 0 b/c reference is invest now
        omega1 = rho/(1-(1/((1+rho)^(N1-0-1))));    %t-1 
        Vt = log(omega0*sum(PVt)/rho);
        Vtp1=log(omega1*(sum(PVtp1)-(1-ejbar*alpha)*eps2)/rho);           %updated in v7 to include learning from others
        deltaLnVj=Vt-Vtp1;                          %   via investment cost uncertainty
        
            
       
        %{
        if mc < 3
            deltaLnVj
            ejbar
        end
        %}

        DELTALnVj(mc,1) = deltaLnVj;
        expmuv(mc,1)=sum(DELTALnVj)/mc;
        j = j + 1;
        end
        
        if errori > limitij | errorj > limitij
            display('Error: sum(PVt) < K or sum(PVtp1) < K more than 100 times.')
            errors = 1;
            j = j + 1;
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

count1=timeseries(expmuv(:,1),1:MC);
count2=timeseries(expsig2v(:,1),1:MC);
plot(count1, '-b')
hold on
plot(count2, '-m')
title('Monte Carlo Simulation')
xlabel('Monte Carlo Iterations')
legend('mean','variance')

end


