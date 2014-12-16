path = fileparts(mfilename('fullpath'))

bs = 15000;

bootstrap = round(100*unifrnd(.01,1,bs,1));
Rdist = zeros(bs,1);
for i = 1: bs
     Rdist(i) = R(bootstrap(i));
end
g = 1

if g == 1;
hist(Rdist,30)
ylabel('Draw Count')
title('Revenue Distribution, 15,000 bootstraps')
xlabel('$ Net Revenue/ Year')

else
figure();
subplot(2,1,1);
hist(Rdist,50)
ylabel('Draw Count')
title('Revenue Distribution, 15,000 bootstraps')

subplot(2,1,2);
hist(R,50)
title('Revenue Distribution, 100 bootstraps')
xlabel('$ Net Revenue/ Year')
ylabel('Draw Count')
hold off

end




