function prob = pricingLewisProb(kappa,deltaT,paramsNMVM,flag)
% This function computes S_T < K at 1year taking advantage of a result of
% the Lewis formula, found in "A SIMPLE OPTION FORMULA FOR GENERAL JUMP-DIFFUSION
% AND OTHER EXPONENTIAL LÉVY PROCESSES" by Alan L. Lewis
%
% INPUT:
% kappa:        the KAPPA of the formula of the paper, not of NIG model
% deltaT:       delta of the function time       
% paramsNMVM:   parameters of the NIG model
% flag:         which model I'm using: 1 = NIG, 2 = VG
%
% OUTPUT
% prob:         wanted probability P(S_1 < K)


% Acquiring parameters of NMVM model
eta = paramsNMVM.eta;
k = paramsNMVM.k;
sigma = paramsNMVM.sigma;
alpha = paramsNMVM.alpha;

if flag == 1 % NIG model
    % Laplace exponent of tempered stable positive random variables
    laplaceExp = @(w) deltaT/k * (1-alpha)/alpha * (1 - (1 + (w.*k.*sigma^2)./(1-alpha)).^alpha); 

elseif flag==2 % VG model
    % Laplace exponent of tempered stable positive random variables
    laplaceExp = @(w) - deltaT/k * log(1 + w.*k.*sigma^2);

end

% Characteristic function of Normal Mean-Variance Mixture model
charFunc = @(z) exp(-1i*z*laplaceExp(eta)).*exp(laplaceExp((z.^2 + 1i*(1 + 2*eta).*z)/2));

% Create the integrand of the function
integrand = @(z) real(exp(1i*z*kappa) .* charFunc(z) ./ (1i*z));

% Compute numerically the integral
int = integral(integrand,0,inf,'AbsTol',1e-13,'RelTol',1e-11);

% Compute probability
prob = 1 - (1/2 + 1/pi * int);

end % pricingLewisProb