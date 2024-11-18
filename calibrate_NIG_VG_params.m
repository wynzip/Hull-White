function [paramsNIG,paramsVG] = calibrate_NIG_VG_params(inputStruct,mktCall,x)
% This function calibrates parameters of the model NIG and VG model
%   
% INPUT:
% inputStruct: contains deltaT, B, F0
% mktCall: call market values
% x: log-moneyness grid
%
% OUTPUT:
% paramsNIG: parameters NIG model (eta,k,sigma)
% paramsVG: parameters VG model (eta,k,sigma)
%
% USES:
% parameters FFT, auxRunPriceF
%


    % Computing the parameters for the FFT grid
    M = 16; % then, N = 2^M
    dz = 0.075; % step of the grid of integral truncation 
    paramsFFT = parametersFFT(M,dz);
    
    % Find the parameters eta,k,sigma of NIG model miniziming the distance
    % between model prices and market prices: MODEL CALIBRATION
    options = optimset('Display', 'Off');
    % eta, kappa, sigma
    inputStruct.alpha = 1/2; % NMVM parameter from text
    calibratedParams = lsqnonlin(@(params) (auxRunPriceF(params(1),params(2),params(3),inputStruct,paramsFFT,x) - mktCall), ...
        [0.001 0.001 0.001], [0 0 0], [1e3 1e3 1e3], [], [], [], [],@(params) nonLinCon(params, alpha), options);
    
    paramsNIG.eta = calibratedParams(1);
    paramsNIG.k = calibratedParams(2);
    paramsNIG.sigma = calibratedParams(3);
    paramsNIG.alpha = 1/2;
    
    % Find the parameters eta,k,sigma of VG model miniziming the distance
    % between model prices and market prices: MODEL CALIBRATION
    options = optimset('Display', 'Off');
    % eta, kappa, sigma
    inputStruct.alpha = 0;
    calibratedParams = lsqnonlin(@(params) (auxRunPriceF(params(1),params(2),params(3),inputStruct,paramsFFT,x) - mktCall), ...
        [0.001 0.001 0.001], [0 0 0], [1e3 1e3 1e3], [], [], [], [],@(params) nonLinCon(params, alpha), options);
    
    % Save the parameters
    paramsVG.eta = calibratedParams(1);
    paramsVG.k = calibratedParams(2);
    paramsVG.sigma = calibratedParams(3);
    paramsVG.alpha = 0;

end % calibrate_NIG_VG_params