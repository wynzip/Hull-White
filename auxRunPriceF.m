function callLewis = auxRunPriceF(eta,k,sigma,inputStruct,paramsFFT,x)
% This is an auxiliary function needed to use the runPricingFourier method
% in the case of exercise 4, where we want to use the FFT pricing to
% calibrate the NMVM, so the parameters eta, k and sigma are unkwowns or
% given as inputs; it's to reconcile code structures of the two exercises
%
% INPUT
% eta           parameter of NMVM model
% k             parameter of NMVM model
% sigma         parameter of NMVM model
% inputStruct   struct containining:
% B             relevant discount factor        
% F0            forward value in t0
% alpha         parameter of NMVM model
% deltaT        parameter of NMVM model
% paramsFFT     struct containing parameters of FFT
% x             moneyness vector
%
% OUTPUT
% callLewis     value of Call using Lewis on moneyness x
%
% USES
% function:     runPricingFourier(B,F0,paramsNMVM,paramsFFT,x,1,1)
%
    
    % Prepare the inputs to send to the actual runPricingFourier function
    B = inputStruct.B;
    F0 = inputStruct.F0;
    paramsNMVM.alpha = inputStruct.alpha;
    paramsNMVM.deltaT = inputStruct.deltaT;

    % Add to the parameters of NMVM model the variables we need to find
    paramsNMVM.eta = eta; 
    paramsNMVM.k = k;
    paramsNMVM.sigma = sigma;

    % Obtain the Call values (flag = 1 for FFT method and 
    % calibrationFlag = 1 because we're indeed calibrating)
    callLewis = runPricingFourier(B,F0,paramsNMVM,paramsFFT,x,1,1);

    
end % auxRunPriceF