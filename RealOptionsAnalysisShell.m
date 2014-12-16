%Begun 11/11/14 by Gregory M. Astill
%This code combines the many modules of code that exist in the folder one
%level above where this file is saved. 
%
%It performs an ex-ante real options analysis of investment values using
%monte-carlo simulation.
%
%Important components: A. Capital Cost distribution, B. Net Revenue
%distribution, C. Neighbor distribution, and D. Environmental Valuation
%distribution.

%A. Capital Cost Distribution
%B. Net Revenue Distribution
%C. Neighbor Distribution
%From KSimulator.m

format compact

datafile = 'C:\Users\Greg\Google Drive\ROA&LearningSpillovers\CowleyData\CowleyDataOut.csv';
data = csvread(datafile,1);

%fid = fopen(datafile);
%VarNames = textscan(fid,'%s',1,'HeaderLines',0)

%Variable Codebook by column
%1  quantiles
%2  pay_personal	
%3  pay_grant	
%4  pay_other	
%5  capau	
%6  anncostau	
%7  ad_perc_income	
%8  methanerevau	
%9  kwcapau	
%10  sell_elect	
%11 elect_wholesalep	
%12 elect_retailp	
%13 coprodau	
%14 instateneigh	
%15 nxtstateneigh

%D. Environmental Valuation Distribution
datafile = 'C:\Users\Greg\Google Drive\ROA&LearningSpillovers\CowleyData\CowleyDataOut_HeterogeneityOfEnvironmentalValuation.csv';
data2 = csvread(datafile,1);

%1.	adopt_sell
%2.	adopt_odor
%3.	adopt_use
%4.	adopt_pr
%5.	adopt_env
%6.	adopt_reg
%7.	adopt_reben
%8.	adopt_ghgben
%9.	adopt_redghg

BS= 1;
N = 1000;
rQuantiles = numel(data(:,1));
cVariables = numel(data(1,:)); 
cVariables2= numel(data2(1,:));
bsK=zeros(BS,2);

for bs=1:BS
            %draws from raw variables using uniform dist
            bsU = zeros(N,cVariables+cVariables2);
            u = unifrnd(0,1,N,cVariables);
            
            %Iterate through N firms
            for U=1:N
                %Iterate through V variables
                for V=1:cVariables
                    %Iterate through Q quantiles
                    for Q=1:rQuantiles
                        if u(U,V) < data(Q,1) 
                            %display yes
                            %u(U,V)
                            %data(Q,1)
                            %Q
                            bsU(U,V) = data(Q,V);
                            break;
                        end
                    end
                end
                %Iterate through V variables of second dataset
                %Rating reasons for AD adoption given by mean and std dev
                %instead of quantiles.
                for V = 1:cVariables2
                            bsU(U,cVariables+ V) = normrnd(data2(2,V),data2(3,V));
                    %pause
                end       
            end
            
            %hist(bsU(:,3));
            %mean = sum(bsU)/N;
            %var   = sum((mean - bsU).^2)/N;
            %change to matrix form
            
            %bsK(bs,:) = [mean var];
            
            
            
end

%hist(bsU(:,cVariables+5))
%pause
%[Kpdf,ri] = ksdensity(bsK(:,1));
%plot(ri, Kpdf)
%pause

K = (100-bsU(:,3)) ./ 100 .* bsU(:,5);
%hist(K)
%6  anncostau	
%7  ad_perc_income	
%8  methanerevau	
%9  kwcapau	
%10  sell_elect	
%11 elect_wholesalep	
%12 elect_retailp	
%13 coprodau
E = 365*24*.8 .* (bsU(:,10) .* bsU(:,11) + (1 - bsU(:,10)) .* bsU(:,12));
C = -bsU(:,6) + bsU(:,8) + E + bsU(:,13);
%14 instateneigh	
%15 nxtstateneigh
NinState = bsU(:,14);
NnextState = bsU(:,15);


%1.	adopt_sell
%5.	adopt_env
EnvValue = bsU(:,cVariables + 1) ./ bsU(:,cVariables + 5) .* C;
%hist(EnvValue)

Kout = [K,C,NinState,NnextState,EnvValue];
savefile = 'Ksim.mat';
save(savefile, 'Kout');






