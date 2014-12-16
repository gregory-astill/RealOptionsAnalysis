%function s = mc_sigma_v10(X, rho, ejbar, alpha, N1, K, Kout, N, MC, graphs) 
    %v5 begun 5/24

%%----------------------------------------------%%       
%%Begin mc_sigma function-----------------------%%
%%----------------------------------------------%%
REnvR, rho, ejbar, alpha, N1, K, Kout, N, MC, graphs

    Np1=N1+1;               %Offset investment lifetime


    PVt = zeros(N1,1);      %Present value vector
    PVtp1= zeros(N1,1);     %Present value of one year delayed investment vector
    muv = zeros(MC,1);      %Mean of present values
    expmuv= zeros(MC,1);    %Expected mean
    mcdeltaLnVjL = 0;       %deltaLNVjL to be saved for each Monte Carlo
    mcdeltaLnVj = zeros(MC,1);
                            %deltaLnVj to be saved for each Monte Carlo
    sig2v = zeros(MC,1);    %variance
    mcsig2vL = 0;           %sig2vL to be saved for each Monte Carlo
    mcsig2v = zeros(MC,1);  %sig2v to be saved for each Monte Carlo
    expsig2v = zeros(MC,1); %Expected variance

    DELTALnVj = zeros(MC,1);%Difference of logged present values.
    errors = 0;             
    mc = 1;
    blah = 0;

        while mc <= MC
            errori = 0;
            errorj = 0;
            limitij= 100;

            TX = zeros(Np1,1);
            tX = zeros(N1,1);
            tp1X=zeros(N1,1);
            u = round(100*unifrnd(0.01,N/100,Np1,1));

                for n1=1:Np1   

                    TX(n1,1) = X(u(n1,1));  
                        %Draws for investment lifetime plus one year

                    if n1 < Np1
                        %Create present value vector for investing today
                        %(Vector elements weighted but not summed)
                        tX(n1,1) = TX(n1,1);
                        PVt(n1,1) = tX(n1,1)/((1+rho)^(n1));
                    end

                    if n1 > 1
                        %Create present value vector for investing tomorrow
                        %(Vector elements weighted but not summed)
                        tp1X(n1-1,1) = TX(n1,1);
                        PVtp1(n1-1,1)=tp1X(n1-1,1)/((1+rho)^(n1-1));
                    end

                end

                omega0 = rho/(1-(1/((1+rho)^(N1-0))));      %t = 0 b/c reference is invest now
                omega1 = rho/(1-(1/((1+rho)^(N1-0-1))));    %t-1 

                Vt  = omega0*sum(PVt)/rho;
                Vtp1= omega1*sum(PVtp1)/rho;
                if Vt <= 0                            %Removed Errors if present value is less than capital cost
                    errori = errori + 1;                    %TODO errors should be present if log(x) below is given x<0
                    display('Error i')
                elseif Vtp1 <= 0
                    errorj = errorj + 1;
                    display('Error j')
                else                    
                    Vt = log(omega0*sum(PVt)/rho);              %Value of investment in perpetuity
                    Vtp1=log(omega1*sum(PVtp1)/rho);            %Value of investment next year in perpetuity
                    deltaLnVj=Vt-Vtp1;                          %Difference of logged investment values
                    eps2 = normrnd(0,abs(deltaLnVj)/3);         %Additional uncertainty that is reduced by learning from others
                    Vnoise = (1-ejbar*alpha)*eps2;              %Learning from others term that goes to zero as more firms adopt
                    deltaLnVj = deltaLnVj + Vnoise;             %Adjusted logged difference in investment values

                     if imag(Vt) ~= 0 | imag(Vtp1) ~= 0
                         warning('Imaginary numbers')
                         pause
                     end

                    DELTALnVj(mc,1) = deltaLnVj;                %Vector of each logged diffence in investment values by iteration
                    expmuv(mc,1)=sum(DELTALnVj)/mc;             %Vector of expected means of logged difference in investment values
                end

                if errori > limitij | errorj > limitij
                    display('Error sum(PVt) < K or sum(PVtp1) < K more than 100 times.')
                        %If too many errors, try again.
                else
                    mc = mc + 1; 
                end
        end
        
        for mc=1:MC

            sig2v(mc,1) = (DELTALnVj(mc,1) - expmuv(MC,1))^2;
            expsig2v(mc,1)=sum(sig2v)/mc;

        end

        expmu = expmuv(MC,1);
        s = expsig2v(MC,1);

        if graphs == 1

            count1=timeseries(expmuv(:,1),1:MC);
            count2=timeseries(expsig2v(:,1),1:MC);
            plot(count1, '-b')
            hold on
            plot(count2, '-m')
            title('Monte Carlo Simulation')
            xlabel('Monte Carlo Iterations')
            legend('mean','variance')

        end


%%----------------------------------------------%%
%%--End mc_sigma function-----------------------%%
%%----------------------------------------------%%



