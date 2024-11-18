function [prob_down_1y,IC_prob_down_1y] = simul_MC_NIG_2y(paramsNMVM,F0,K,nSim,deltaT)
% This function checks evolution of the underlying NIG in 1 year and
% computes probability of having S_T < K at 1 year
%
%   INPUT
%   F0:             fwd value in 0
%   nSim:           number of MC simulations
%   deltaT:     needed yearfrac
%   paramsNMVM:     parameters of the NIG model
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
    
    % Select the indexes of the ones below the threshold
    idx_down = forwardValueT < K;
    
    % Compute the probability at 1y and its confidence interval (normfit)
    [prob_down_1y,~,IC_prob_down_1y] = normfit(idx_down);
    

    figure()
    histogram(forwardValueT, 'FaceColor', "#4DBEEE", 'EdgeColor', 	"#0072BD")
    hold on
    xline(K, 'LineWidth', 2, 'Color', "#A2142F")
    title('Forward Probability Distribution - MC')
    hold off
    

end % simul_MC_NIG_2y







