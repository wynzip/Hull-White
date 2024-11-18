function [I] = integralQuad(charFunc,x,absTol,relTol)

% This function computes the integral of the Lewis formula using the
% Quadrature method (numerical approximation of the integral)
% INPUT
% charFunc      characteristic function of the model used for the underlying
% x             moneyness vector on which to evaluate integral
% absTol        absolute tolerance as numerical integration parameter
% relTol        relative tolerance as numerical integration parameter
%
% OUTPUT
% I             values of the integral computed
%


    % Create function handle of integrand of Lewis integral (missing the
    % exponential term)
    funLewisIncomplete = @(z) 1/(2*pi) .* 1./(z.^2 + 1/4) .* charFunc(-z - 1i/2);
    
    % Preallocate vector of integral values (one for each moneyness value)
    intQuad = zeros(length(x),1);
    
    for n=1:length(x) % loop over all the moneyness values

        currMoneyness = x(n); % current moneyness value
        
        % Now we can write the entire integrand of the Lewis formula, as it
        % depends from the chosen moneyness value
        funLewis = @(z) funLewisIncomplete(z) .* exp(-1i*currMoneyness.*z);

        % Numerical evaluation of the Lewis formula integral
        intQuad(n) = integral(funLewis,-inf,inf,'AbsTol',absTol,'RelTol',relTol); % SET PRECISION OF INTEGRAL
    end

     % Now let's check that the imaginary part is negligible, as wanted
    if max(abs(imag(intQuad))) < absTol
        I = real(intQuad);
    else
        I = zeros(length(xGrid),1); % we put zeros so we know the result is wrong
        disp('Problem in Quadrature computation: imaginary part is NON-negligible')
    end


end % integralQuad
