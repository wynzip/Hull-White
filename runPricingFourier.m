function callLewis = runPricingFourier(B,F0,paramsNMVM,paramsFFT,x,flag,calibFlag,absTol,relTol)
   
% This is a function to price EU Calls on the moneyness values x using the 
% Lewis formula and two possible methods, FFT or Quadrature
%
% INPUT
% B             relevant discount factor        
% F0            forward value in t0
% paramsNMVM    struct containing parameters of NMVM model
% paramsFFT     struct containing parameters of FFT
% x             moneyness vector
% flag          method: 1 = FFT, 2 = Quadrature
% calibFlag     1 = calibration is ongoing, 0 = otherwise
% absTol        absolute tolerance to use in Quad method
% relTol        relative tolerance to use in Quad method
%
% OUTPUT
% callLewis     value of Call using Lewis on moneyness x
%
% USES
% function      integralFFT(paramsFFT,charFunc,calibFlag)
%               integralQuad(charFunc,x,absTol,relTol)
%
    
    % Acquiring parameters of NMVM model
    alpha = paramsNMVM.alpha;
    sigma = paramsNMVM.sigma;
    k = paramsNMVM.k;
    eta = paramsNMVM.eta;
    deltaT = paramsNMVM.deltaT;
    if alpha == 0
          % Laplace exponent of tempered stable positive random variables
          laplaceExp = @(w) - deltaT/k * log(1 + w.*k.*sigma^2); 
    else
          % Laplace exponent of tempered stable positive random variables
          laplaceExp = @(w) deltaT/k * (1-alpha)/alpha * (1 - (1 + (w.*k.*sigma^2)./(1-alpha)).^alpha); 
    end

    % Characteristic function of Normal Mean-Variance Mixture model
    charFunc = @(z) exp(-1i*z*laplaceExp(eta)).*exp(laplaceExp((z.^2 + 1i*(1 + 2*eta).*z)/2));

    switch flag % method: 1 = FFT, 2 = Quadrature

        case 1
            % Compute the integral in the Lewis formula using the FFT method   
            [xGrid,I] = integralFFT(paramsFFT,charFunc,calibFlag);

            % Interpolate the integral on the wanted moneyness values
            intRequest = interp1(xGrid,I,x);
            
        case 2
            % Compute the integral in the Lewis formula using the Quadrature method
            intRequest = integralQuad(charFunc,x,absTol,relTol); 
            % already computed on the wanted moneyness grid

    end

    % Now we can price the Call option on the wanted moneyness values
    callLewis = B*F0*(1 - exp(-x/2).*intRequest);


end % runPricingFourier