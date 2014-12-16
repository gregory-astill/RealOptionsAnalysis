function trigger = triggersolve(rho, N1, K)

trigger_sum = zeros(1,N1);

    for i=1:N1
        trigger_sum(1,i) = 1/(1+rho)^i;
    end

trigger = fzero(@(x) K-sum(trigger_sum*x), 50);

end