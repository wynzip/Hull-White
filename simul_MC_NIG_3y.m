function [prob_down_down,prob_down_up,prob_up_down,prob_up_up] = simul_MC_NIG_3y(paramsNMVM,F0,K,nSim,deltaT_vec)
% This function checks evolution of the underlying NIG in 2 years and
% computes probabilities for each of the 4 scenarios of the product
%
%   INPUT
%   F0:             fwd value in 0
%   nSim:           number of MC simulations
%   deltaT_vec:     vec of needed yearfracs
%   paramsNMVM:     parameters of the NIG model
%
%   OUTPUT
%   prob_down_down: probability of having S_1 < K and S_2 < K
%   prob_down_up:   probability of having S_1 < K and S_2 > K
%   prob_up_down:   probability of having S_1 > K and S_2 < K
%   prob_up_up:     probability of having S_1 > K and S_2 > K
%   
    
    % Set random seed for reproducibility of results
    rng(51)
    
    % Acquiring parameters of NMVM model
    eta = paramsNMVM.eta;
    k = paramsNMVM.k;
    sigma = paramsNMVM.sigma;
    
    % Extract deltaT first year
    deltaT = deltaT_vec(1);
    
    % Set up parameters of our InverseGaussian distribution
    muIG = 1; % scale parameter of IG
    lambdaIG = deltaT/k; % shape parameter of IG
    
    % Simulate random values following the Inverse Gaussian distribution
    valuesIG = random('InverseGaussian',muIG,lambdaIG,[nSim,1]);
    % Simulate random values following gaussian distribution
    valuesGauss = randn(nSim,1);
    
    % Laplace exponent of tempered stable positive random variables
    laplaceExpNIG = @(w) deltaT/k * (1 - sqrt(1 + 2*w.*k.*sigma^2));
    % Compute the value of mu
    mu = - laplaceExpNIG(eta); 
    
    % Simulate ft value for the 1st year following the model
    ft = mu - sigma^2*(0.5 + eta)*deltaT*valuesIG + sigma*sqrt(deltaT)*sqrt(valuesIG).*valuesGauss;
    
    % Here follows F(t,t) = F0*exp(ft), the forward price
    forwardValueT = F0*exp(ft);
    
    % Select indexes of the ones above/below the threshold at 1y (we need them later)
    idx_up_1y = forwardValueT > K;
    idx_down_1y = forwardValueT < K;
    
    % Extract deltaT
    deltaT = deltaT_vec(2) - deltaT_vec(1);
    
    % Set up parameters of our InverseGaussian distribution
    muIG = 1; % scale parameter of IG
    lambdaIG = deltaT/k; % shape parameter of IG
    
    % Simulate random values following the Inverse Gaussian distribution
    % we want to have Nsim simulations for every moneyness grid value
    valuesIG = random('InverseGaussian',muIG,lambdaIG,[nSim,1]);
    
    % Simulate random values following gaussian distribution
    valuesGauss = randn(nSim,1);
    
    % Laplace exponent of tempered stable positive random variables
    laplaceExpNIG = @(w) deltaT/k * (1 - sqrt(1 + 2*w.*k.*sigma^2));
    
    % Compute the value of mu
    mu = - laplaceExpNIG(eta); 
    
    % Simulate ft value following the model
    ft = ft + mu - sigma^2*(0.5 + eta)*deltaT*valuesIG + sigma*sqrt(deltaT)*sqrt(valuesIG).*valuesGauss;
    
    % Here follows F(t,t) = F0*exp(ft), the forward price
    forwardValueT = F0*exp(ft);
    
    % Check the underlying situation at 2 years:
    idx_up_2y = forwardValueT > K;
    idx_down_2y = forwardValueT < K;

    % Obtain final four scenarios
    idx_down_down = idx_down_1y.*idx_down_2y;
    idx_down_up = idx_down_1y.*idx_up_2y;
    idx_up_down = idx_up_1y.*idx_down_2y;
    idx_up_up = idx_up_1y.*idx_up_2y;
    
    % Obtain the final 4 needed probabilities
    prob_down_down = sum(idx_down_down)/nSim;
    prob_down_up = sum(idx_down_up)/nSim;
    prob_up_down = sum(idx_up_down)/nSim;
    prob_up_up = sum(idx_up_up)/nSim;
    
end % simul_MC_NIG_3y







