function [xGrid,I] = integralFFT(paramsFFT,charFunc,calibFlag)
  
% This function computes the integral of the Lewis formula using the FFT
% method
% INPUT
% paramsFFT     struct containing parameters of FFT
% charFunc      characteristic function of the model used for the underlying
% calibFlag     1 = calibration is ongoing, 0 = otherwise
%
% OUTPUT
% xGrid         moneyness grid on which integral is computed
% I             values of the integral computed
%

    % Exctracting from the struct all the needed parameters
    zGrid = paramsFFT.zGrid; % integration grid (truncation)
    dz = paramsFFT.dz; % step of the integration grid
    x1 = paramsFFT.x1; % lowest point of moneyness grid
    xGrid = paramsFFT.xGrid; % moneyness grid
    N = paramsFFT.N; % number of points in each grid
    z1 = paramsFFT.z1; % lowest point of integration grid (truncation)
    
    % Create function handle of f(z) of Lewis integral --> for Fourier transform
    funLewis = @(z) 1/(2*pi) .* 1./(z.^2 + 1/4) .* charFunc(-z -1i/2); 
    funLewisGrid = funLewis(zGrid); % valued f(Zj) for every Zj (grid of integration)
    
    % We have f(Zj), we have to apply the phase to it, for the FFT method
    j = (0:N-1)';
    funFFT = funLewisGrid .* exp(-1i*x1*dz.*j); % this is the function to be sent inside the FFT algo
    
    % Computing the prefactor we will need to multiply to the FFT result
    prefactor = dz.*exp(-1i*z1.*xGrid); % one for each moneyness value
    
    % Compute integral using the FFT method, one for each moneyness value
    intFFTGrid = prefactor.*fft(funFFT); 

    % Now let's check that the imaginary part is negligible, as wanted
    if calibFlag == 0  % NO calibration
        if max(abs(imag(intFFTGrid))) < 1e-10 
            I = real(intFFTGrid);
        else
            I = zeros(length(xGrid),1); % we put zeros so we know the result is wrong
            disp('Problem in FFT computation: imaginary part is NON-negligible')
        end
    else % calibration
        % The rationale is: if we're using this function for the Model
        % Calibration, especially for the first values of the parameters
        % we'll have integrals with rather big imaginary parts, because
        % we're using wrong parameters, so we don't check the imaginary
        % part in the model calibration case
        I = real(intFFTGrid);
    end

end % integralFFT