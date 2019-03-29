function sample = Metropolis_Hastings(x0, log_density, random_step, N_sample, noisy)
N = length(x0); sample = zeros(N, N_sample);
sample(:,1) = x0; x = x0;
l_x = log_density(x); k = 1; ks = 1;
while ks < N_sample
    w = normrnd(0, 1, N, 1);
    y = x + random_step*w;
    l_y = log_density(y);
    
    if noisy
        l_x = log_density(x);
    end
    
    ks = ks+1;
    if l_y - l_x > log(rand)
        x = y;
        l_x = l_y;
        k = k+1;
        sample(:,ks) = x;
    else
        sample(:,ks) = x;
    end
    accrate = k/ks*100;
    
    if mod(ks, 10) == 0
        disp(['iteration        ', num2str(ks)])
        disp(['acceptance ratio ', num2str(accrate)]);
        disp(['current estimate ', num2str(sample(:, ks)')])
    end
end