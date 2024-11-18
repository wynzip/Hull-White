function [x_tree,stoc_disc_mat,yearly_fwd_disc,dt,tree_yearfracs,yf_365_yearly] = initialize_Tree(M,dates,discounts,a,sigma)
% This function initializes the trinomial HW tree, computing some needed
% variables
%
% INPUT:
% M:                number of time steps in the tree (from 0 to 10 years)
% dates:            dates of bootstrap
% discounts:        discount factors of bootstrap
% a:                parameter of HW model
% sigma:            parameter of HW model
%
% OUTPUT:
% x_tree:            values of the x_i nodes on the tree
% stoc_disc_mat:     "tree"(matrix) of stochastic discounts D(0,t_i,t_i+dt)
% yearly_fwd_disc:   matrix of yearly fwd disc B(0,year,following_years)
% dt:                time step in the tree
% tree_yearfracs:    time grid of the tree
% yf_365_yearly:     yearly yearfracs in Act/365
%

    C_ACT_365 = 3;
    %% Tree H&W: let's start setting up the parameters
    
    % End date of the contract
    date_end = datenum('2018/02/19');
    yf_end_365 = yearfrac(dates(1),date_end,C_ACT_365);
    
    % Find the time grid of the tree
    tree_yearfracs = linspace(0,yf_end_365,M+1)'; % time grid of the tree
    
    % Find the known forward discounts B(0,ti,ti+1) in t0 using the bootstrap
    fwd_discounts = disc_factor_fwd(dates, discounts, tree_yearfracs(2:end));
    
    % Parameters HW and tree
    dt = tree_yearfracs(2) - tree_yearfracs(1); % time space
    mu_hat = 1 - exp(-a*dt);
    
    % Bounds of the tree
    l_max = (1-sqrt(2/3))/mu_hat;
    l_max = ceil(l_max); % l_min = -l_max;
    
    % Other params
    sigma_hat = sigma*sqrt((1-exp(-2*a*dt))/(2*a));
    dx = sqrt(3)*sigma_hat;
    sigma_hat_star = sigma/a*sqrt(dt - 2*(1 - exp(-a*dt))/a + (1 - exp(-2*a*dt))/(2*a));
    
    %% Creating the matrix of X dynamics
    x0 = 0;
    x_tree = zeros(l_max*2+1, M+1);
    x_tree(l_max+1, :) = x0;
    % x_final = x0 + delta_x*(l_min:l_max)';
    
    for i = 1:l_max
        % Down-values
        x_tree(l_max+1+i, (i+1):end) = -i*dx*ones(M+1-i, 1)';
        % Up-values
        x_tree(l_max+1-i, (i+1):end) = +i*dx*ones(M+1-i, 1)';
    end
    
    %% Creating the matrix of simulated spot discount factors in the tree, from Lemma 1 corollary
    
    % Define the function handle of the HJM volatility sigma(t1,t2)
    sigma_HJM = @(t1,t2) sigma/a*(1 - exp(- a*(t2 - t1)));
    
    % Evaluate HJM sigma(0,dt)
    sigma_0_dt = sigma_HJM(0,dt);
    
    % Find the integrand part of the formula
    integrand = @(u,t) sigma_HJM(u,t+dt).^2 - sigma_HJM(u,t).^2;
    
    % Initialize the integrals vector
    int_vec = zeros(length(tree_yearfracs),1);
    
    % Evaluate the integral for each ti - time step on the tree
    for i=1:length(tree_yearfracs)
    
        integrand_bis = @(u) integrand(u,tree_yearfracs(i));
    
        % Evaluate each integral for ti
        int_vec(i) = integral(integrand_bis,0,tree_yearfracs(i));
    end
    
    % Initialize spot discount matrix
    spot_disc_mat = zeros(size(x_tree));
    
    % Fill the discount matrix for the easy values
    for i = length(tree_yearfracs)-1:-1:l_max+1
    
        % Compute the needed exponential term of the formula
        exp_term = exp(-x_tree(:,i)*sigma_0_dt/sigma - 1/2*int_vec(i));
        
        % Compute the discount B(ti,ti,ti+dt)
        spot_disc_mat(:,i) = fwd_discounts(i)*exp_term;
    
    end
    
    for i = l_max:-1:1
    
        % Find the non-zero values of x in the tree
        x = x_tree(l_max+1-i+1:l_max+1+i-1,i);
        % the indexes inside are a bit messy but they work, they need to be
        % symmetric around the l_max+1 value which is the middle tree value
    
        % Compute the needed exponential term of the formula
        exp_term = exp(-x*sigma_0_dt/sigma - 1/2*int_vec(i));
    
        % Find and update the new discount factor
        spot_disc_mat(l_max+1-i+1:l_max+1+i-1,i) = fwd_discounts(i)*exp_term;
    
    end
    
    %% Creating the matrix of simulated stochastic discount factors in the tree, from Lemma 1
    
    % Initialize stochastic discounts matrix
    stoc_disc_mat = zeros(size(x_tree));

    % initialize the exponential term
    exp_term = zeros(size(x_tree,1),3);

    %compute probabilities
    [p_up, p_m, p_down] = compute_probabilities_tree(l_max, mu_hat);
    
    % Fill the discount matrix for the easy values
    for i = length(tree_yearfracs):-1:l_max+1
        
        % compute exponential term for the l_max scenario
        exp_term(1,1) = exp(-0.5*sigma_hat_star^2 - sigma_hat_star/sigma*( mu_hat*x_tree(1,i)));
        exp_term(1,2) = exp(-0.5*sigma_hat_star^2 - sigma_hat_star/sigma*(-dx+ mu_hat*x_tree(1,i)));
        exp_term(1,3) = exp(-0.5*sigma_hat_star^2 - sigma_hat_star/sigma*(-2*dx+mu_hat*x_tree(1,i)));
        % compute exponential term for the l scenario
        exp_term(2:end-1,1) = exp(-0.5*sigma_hat_star^2 - sigma_hat_star/sigma*( dx +mu_hat*x_tree(1,i)));
        exp_term(2:end-1,2) = exp(-0.5*sigma_hat_star^2 - sigma_hat_star/sigma*( mu_hat*x_tree(1,i)));
        exp_term(2:end-1,3) = exp(-0.5*sigma_hat_star^2 - sigma_hat_star/sigma*(-dx +mu_hat*x_tree(1,i)));
        % compute exponential term for the l_min scenario
        exp_term(end,1) = exp(-0.5*sigma_hat_star^2 - sigma_hat_star/sigma*( +2*dx +mu_hat*x_tree(1,i)));
        exp_term(end,2) = exp(-0.5*sigma_hat_star^2 - sigma_hat_star/sigma*( dx+mu_hat*x_tree(1,i)));
        exp_term(end,3) = exp(-0.5*sigma_hat_star^2 - sigma_hat_star/sigma*( mu_hat*x_tree(1,i)));
        % exponential term to be used in the formula
        exp_term = exp_term(:,1).*p_up + exp_term(:,2).*p_m + exp_term(:,3).*p_down;
        
        % Compute the discount B(ti,ti,ti+dt)
        stoc_disc_mat(:,i) = spot_disc_mat(:,i).*exp_term;
    
    end
    
    for i = l_max:-1:1
    
        % Find the non-zero values of x in the tree
        x = x_tree(l_max+1-i+1:l_max+1+i-1,i);

        %initialize exponential term
        exp_term = zeros(length(x), 3);

        %cut the probabilities to the needed values
        p_up_temp = p_up(l_max+1-i+1:l_max+1+i-1);
        p_m_temp = p_m(l_max+1-i+1:l_max+1+i-1);
        p_down_temp = p_down(l_max+1-i+1:l_max+1+i-1);
        % the indexes inside are a bit messy but they work, they need to be
        % symmetric around the l_max+1 value which is the middle tree value
    
        % compute exponential term for the l scenario
        exp_term(:,1) = exp(-0.5*sigma_hat_star^2 - sigma_hat_star/sigma*( dx +mu_hat*x));
        exp_term(:,2) = exp(-0.5*sigma_hat_star^2 - sigma_hat_star/sigma*( mu_hat*x));
        exp_term(:,3) = exp(-0.5*sigma_hat_star^2 - sigma_hat_star/sigma*(-dx +mu_hat*x));
        % exponential term for the formula
        exp_term = exp_term(:,1).*p_up_temp + exp_term(:,2).*p_m_temp + exp_term(:,3).*p_down_temp;

        % Find and update the new discount factor
        stoc_disc_mat(l_max+1-i+1:l_max+1+i-1,i) = spot_disc_mat(l_max+1-i+1:l_max+1+i-1,i).*exp_term;
    
    end
    
    %% Precompute all the needed fwd discount yearly for the tree
    
    % Number of years of contract
    years = 10;
    % yearly yearfracs
    yf_365_yearly = linspace(0,yf_end_365,years+1)';
    yf_365_yearly(1) = [];
    
    yearly_fwd_disc = zeros(years-2,years-2);
    
    % obtaining zero rates from bootstrapped data
    yf = yearfrac(dates(1), dates(2:end), C_ACT_365);  %delta( t_0, t_market)
    zero_rates= -log(discounts(2:end))./yf; %zero rates from market
    
    zero_rates=[zero_rates(1); zero_rates];
    yf=[0; yf];
    
    for i = years-1:-1:2
        
        % temporary needed yearfracs for fwd discount
        year_frac = yf_365_yearly(i:end);
        
        % interpolating on zero rates
        new_zero_rates = interp1(yf, zero_rates, year_frac, 'linear'); % zero rates for payment dates
        
        %discount factors corresponding to dates of payments
        spot_discount = exp(-new_zero_rates .* year_frac); 
    
        % compute fwd discounts
        yearly_fwd_disc(i-1,1:years-i) = (spot_discount(2:end)./spot_discount(1))';
    
    end


end % initialize_Tree