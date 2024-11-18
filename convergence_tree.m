function convergence_tree(dates,discounts)
%plot the prices of the bermudan option varying the number of time steps in
%the tree
% 
% INPUT
% dates:        vector of dates
% discounts:    vector of discounts
%

temp = [25:25:3650];
MM=[ 10, temp];

prices = ones(length(MM),1);

for i=1:length(MM)
    M=MM(i);

    % Parameters of Hull and White model
    a = 0.11;
    sigma = 0.008;
    
    % Initialize needed parameters of the tree
    [x_tree,stoc_disc_mat,yearly_fwd_disc,dt,tree_yearfracs,yf_365_yearly] = initialize_Tree(M,dates,discounts,a,sigma);
    
    % Mu_hat for the OU process
    mu_hat = 1 - exp(-a*dt);
    % Recompute bounds of the tree
    l_max = (1-sqrt(2/3))/mu_hat;
    l_max = ceil(l_max);
    years = 10;
    
    % Compute bermudan swaption price thanks to the tree
    prices(i) = bermudan_tree(x_tree,stoc_disc_mat,yearly_fwd_disc,dt,tree_yearfracs,yf_365_yearly,mu_hat,l_max,years,M,dates,a,sigma);
end

figure()
plot(MM, prices, 'LineWidth', 2)
title('Convergence of the Trinomial Tree')
xlabel('Number of time steps')
ylabel('Price of the Bermudan Option')


end %convergence_tree