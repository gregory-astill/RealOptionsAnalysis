function s = mc_sigma_v11(X, rho, ejbar, alpha, N1, N, MC, graphs, eps) 
    %v5 begun 5/24
    %v10 restructed error catching, removed unnecessary variables
    %v10 removed tX and tp1X
    %Removing for loops to increase code speed.

%%----------------------------------------------%%       
%%Begin mc_sigma function-----------------------%%
%%----------------------------------------------%%
    Np1=N1+1;               %Offset investment lifetime


    PVt = zeros(N1,1);      %Present value vector
    PVtp1= zeros(N1,1);     %Present value of one year delayed investment vector
    expmuv= zeros(MC,1);    %Expected mean
    %mcdeltaLnVjL = 0;       %deltaLNVjL to be saved for each Monte Carlo
    %mcdeltaLnVj = zeros(MC,1);
                            %deltaLnVj to be saved for each Monte Carlo
                            
    Vt = zeros(MC,1);
    Vtp1 = zeros(MC,1);
    sig2v = zeros(MC,1);    %variance
    %mcsig2vL = 0;           %sig2vL to be saved for each Monte Carlo
    %mcsig2v = zeros(MC,1);  %sig2v to be saved for each Monte Carlo
    expsig2v = zeros(MC,1); %Expected variance

    DELTALnVj = zeros(MC,1);%Difference of logged present values.         
    mc = 1;

        while mc <= MC
            errori = 0;
            errorj = 0;
            limitij= 100;

            u = round(100*unifrnd(0.01,N/100,Np1,1));

            %Draws for investment lifetime plus one year
            TX = X(u);
            denom = (1+rho).^(1:Np1-1);
            %Create present value vector for investing today
            %(Vector elements weighted but not summed)
            %tX(n1,1) = TX(n1,1);
            PVt = TX(1:Np1-1)./denom';
            %Create present value vector for investing tomorrow
            %(Vector elements weighted but not summed)
            %tp1X(n1-1,1) = TX(n1,1);
            PVtp1 = TX(2:Np1)./denom';

                omega0 = rho/(1-(1/((1+rho)^(N1-0))));      %t = 0 b/c reference is invest now
                omega1 = rho/(1-(1/((1+rho)^(N1-0))));    %t-1 

                Vt(mc,1)  = omega0*sum(PVt)/rho;
                Vtp1(mc,1)= omega1*sum(PVtp1)/rho;
                if Vt(mc,1) <= 0                            %Removed Errors if present value is less than capital cost
                    errori = errori + 1;                    %TODO errors should be present if log(x) below is given x<0
                    display('Error i')
                elseif Vtp1(mc,1) <= 0
                    errorj = errorj + 1;
                    display('Error j')
                else                    
                    Vtlog  = log(omega0*sum(PVt)/rho);              %Value of investment in perpetuity
                    Vtp1log= log(omega1*sum(PVtp1)/rho);            %Value of investment next year in perpetuity
                    deltaLnVj=Vtlog-Vtp1log;                          %Difference of logged investment values
                    eps2 = normrnd(0,abs(deltaLnVj)/eps);         %Additional uncertainty that is reduced by learning from others
                    %eps2 = 0;                                   %No noise for now
                    Vnoise = (1-ejbar*alpha)*eps2;              %Learning from others term that goes to zero as more firms adopt
                    deltaLnVj = deltaLnVj + Vnoise;             %Adjusted logged difference in investment values

                    DELTALnVj(mc,1) = deltaLnVj;                %Vector of each logged diffence in investment values by iteration
                    expmuv(mc,1)=sum(DELTALnVj)/mc;             %Vector of expected means of logged difference in investment values
                end

                if errori > limitij || errorj > limitij
                    display('Error sum(PVt) < K or sum(PVtp1) < K more than 100 times.')
                        %If too many errors, try again.
                else
                    mc = mc + 1; 
                end
        end
        
        sig2v = (DELTALnVj - expmuv(size(DELTALnVj,2),1)).^2;    
        denom = [1:size(DELTALnVj,2)];
        expsig2v = cumsum(sig2v)./denom;

        s = expsig2v(MC,1);

        if graphs == 1
            figure('OuterPosition', [0,30,1280,760]);
            
            subplot(2,4,1);
            hold on
            hist(Vt*rho/omega0,MC/50);
            title('Present Value')
            line([mean(Vt*rho/omega0) mean(Vt*rho/omega0)],[0,100], 'LineWidth', 0.9, 'Color', [1 0 1])
            xminxmax =xlim;
            
            subplot(2,4,5);
            hold on
            hist(Vtp1*rho/omega1,MC/50);
            title('Present Value t + 1')
            xlim([xminxmax(1),xminxmax(2)])
            line([mean(Vtp1*rho/omega1) mean(Vtp1*rho/omega1)],[0,100], 'LineWidth', 0.9, 'Color', [1 0 1])
     
            subplot(2,4,2);
            hold on
            hist(Vt,MC/50);
            title('Non-logged Value at t')
            line([mean(Vt) mean(Vt)],[0 100], 'LineWidth', 0.9, 'Color', [1 0 1])
            xminxmax = xlim;
            
            subplot(2,4,6);  
            hold on
            hist(Vtp1,MC/50);
            title('Non-logged Value at t+1')
            xlim([xminxmax(1),xminxmax(2)])
            line([mean(Vtp1) mean(Vtp1)],[0 100], 'LineWidth', 0.9, 'Color', [1 0 1])
            hold off
            
            subplot(2,4,3);
            count1=timeseries(DELTALnVj,1:MC);
            plot(count1, '-b')
            title('Delta Log Diff of Values')
            xlabel('MC Iterations')
            
            subplot(2,4,7);
            hold on
            hist(DELTALnVj,MC/50);
            title('Delta Log Diff of Values')
            line([mean(DELTALnVj) mean(DELTALnVj)],[0 100], 'LineWidth', 0.9, 'Color', [1 0 1])
            
            subplot(2,4,4);           
            count1=timeseries(expmuv(:,1),1:MC);
            plot(count1, '-b')
            title('MC Mean Estimation')
            xlabel('Monte Carlo Iterations')
            legend('mean')
            
            subplot(2,4,8);
            count2=timeseries(expsig2v(:,1),1:MC);          
            plot(count2, '-m')
            title('MC Variance Estimation')
            xlabel('Monte Carlo Iterations')
            legend('variance')
            
        end


%%----------------------------------------------%%
%%--End mc_sigma function-----------------------%%
%%----------------------------------------------%%



