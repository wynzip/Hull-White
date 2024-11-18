function [prob_down_1y,IC_prob_down_1y] = simul_MC_VG(paramsNMVM,F0,K,nSim,deltaT)
% This function checks evolution of the underlying VG in 1 year and
% computes probability of having S_T < K at 1 year
%
%   INPUT
%   F0:             fwd value in 0
%   nSim:           number of MC simulations
%   deltaT:         needed yearfrac
%   paramsNMVM:     parameters of the VG model
%
%   OUTPUT:
%   prob_down_1y:   probability of having S_T < K at 1 year
%  

    % Set random seed for reproducibility of results
    rng(51)

    % Acquiring parameters of NMVM model
    eta = paramsNMVM.eta;
    k = paramsNMVM.k;
    sigma = paramsNMVM.sigma;

    % Set up parameters of our InverseGaussian distribution
    aG = deltaT/k; % scale parameter of IG
    bG = k/deltaT; % shape parameter of IG

    % Simulate random values following the Inverse Gaussian distribution
    % we want to have Nsim simulations for every moneyness grid value
    valuesGamma = random('Gamma',aG,bG,[nSim,1]);
    
    % Simulate random values following gaussian distribution
    valuesGauss = randn(nSim,1);
    
    % Laplace exponent of tempered stable positive random variables
    laplaceExpVG = @(w) - deltaT/k * log(1 + w.*k.*sigma^2);
    
    % Compute the value of mu
    mu = - laplaceExpVG(eta); 

    % Simulate ft value following the model
    ft = mu - sigma^2*(0.5 + eta)*deltaT*valuesGamma + sigma*sqrt(deltaT)*sqrt(valuesGamma).*valuesGauss;

    % Here follows F(t,t) = F0*exp(ft), the forward price
    forwardValueT = F0*exp(ft);
    
    % Compute the probability
    idx_down = forwardValueT < K;

    % Compute the probability at 1y and its confidence interval (normfit)
    [prob_down_1y,~,IC_prob_down_1y] = normfit(idx_down);

end % simul_MC_VG