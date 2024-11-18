function [price_european] = european_tree(x_tree,stoc_disc_mat,yearly_fwd_disc,dt,tree_yearfracs,mu_hat,l_max,M,dates,a,sigma)
% This function computes the european Swaption price using the trinomial 
% HW tree (fwd starting swap 9-10 years)
%
% INPUT:
% x_tree:            values of the x_i nodes on the tree
% stoc_disc_mat:     "tree"(matrix) of stochastic discounts D(0,t_i,t_i+dt)
% yearly_fwd_disc:   matrix of yearly fwd disc B(0,year,following_years)
% dt:                time step in the tree
% tree_yearfracs:    time grid of the tree
% yf_365_yearly:     yearly yearfracs in Act/365
% mu_hat:            mu_hat for the OU process   
% l_max:             max bound of the tree
% years:             number of years of the contract
% M:                 number of time steps in the tree (from 0 to 10 years)
% dates:            dates of bootstrap
% a:                parameter of HW model
% sigma:            parameter of HW model
%
% OUTPUT:
% price_europeann:  price of the european swaption computed
%

%% Let's start the tree with the 9th year payoff evaluation
C_30_360 = 6;
% Strike of the bermudan swaption
K = 0.05;

% Find the time indexes of the years on the tree
idxs_years = round(1:(1/dt):M+1)';
idxs_years(1) = []; % we don't need the index of year 0

% Precompute all the yearly yearfracs for the BPV computation (30/360)
date_end = datenum('2018/02/19');
yf_end_30_360 = yearfrac(dates(1),date_end,C_30_360);
yf_bpv_yearly = linspace(0,yf_end_30_360,M+1)';
yf_bpv_yearly = yf_bpv_yearly(idxs_years);

% Compute the needed spot discounts B(9y,10y) for each node at 9y
curr_year = 9;
idx_year = idxs_years(curr_year);

% Define the function handle of the HJM volatility sigma(t1,t2)
sigma_HJM = @(t1,t2) sigma/a*(1 - exp(- a*(t2 - t1)));

% Find the integrand part of the formula
integrand = @(u,t1,tau) sigma_HJM(u,t1+tau).^2 - sigma_HJM(u,t1).^2;

% Vector of nodes at 9th year
x_nodes_9y = x_tree(:,idx_year);

% Compute the yearly spot discount by steps
yearly_spot_disc_9y = yearly_fwd_disc(end,1)*ones(length(x_nodes_9y),1);

% integral term
t_i = tree_yearfracs(idx_year);
tau = tree_yearfracs(idxs_years(curr_year+1)) - tree_yearfracs(idx_year);
int = integral(@(u) integrand(u,t_i,tau),0,t_i);
% exponential term
exp_term = exp(-x_nodes_9y*sigma_HJM(0,tau)/sigma - 1/2*int);

yearly_spot_disc_9y = yearly_spot_disc_9y.*exp_term; % multiply everything

% Compute the bpv at 9y-10y
bpv_9y_10y = (yf_bpv_yearly(curr_year+1) - yf_bpv_yearly(curr_year))*yearly_spot_disc_9y;

% Compute the swap rate 9y-10y
swap_rate_9y = (1 - yearly_spot_disc_9y)./bpv_9y_10y;

payoff_9y = bpv_9y_10y.*max(swap_rate_9y - K, 0);

%% Let's start the actual tree

curr_year = 9; % we start from year 9

%initialize expected payoff for the tree
expected_payoff = payoff_9y;

% Prepare tree probabilities
[p_up, p_m, p_down] = compute_probabilities_tree(l_max, mu_hat);

for j = idxs_years(curr_year)-1:-1:1 % iterate over time steps

    % Discount the values
    expected_payoff(1) = p_up(1)*expected_payoff(1) + p_m(1)*expected_payoff(2) ...
        + p_down(1)*expected_payoff(3);

    expected_payoff(2:end-1) = p_up(2:end-1).*expected_payoff(1:end-2) + ...
        p_m(2:end-1).*expected_payoff(2:end-1) + p_down(2:end-1).*expected_payoff(3:end) ;

    expected_payoff(end) = p_up(end)*expected_payoff(end-2) + ...
        p_m(end)*expected_payoff(end-1) + p_down(end)*expected_payoff(end);

    expected_payoff = stoc_disc_mat(:, j) .* expected_payoff;
   
end

% Select the computed price of the bermudan swaption in t0
price_european = expected_payoff(l_max+1);

end