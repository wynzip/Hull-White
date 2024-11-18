function swaptionPrices = swaptionPriceJamshidian (jData)
% This function computes the swaption prices from 2y-9y to 10y
% via Jamshidian approach 
%
% INPUT:
% jData:            struct with:
%                    - a:                Hull-White a
%                    - sigma:            Hull-White sigma
%                    - strike:           coupon bond put strike
%                    - delta_YF365:      delta from t0 to payment dates 
%                                         with ACT 365
%                    - spot_discount:    spot discount factor vector   
%                    - c_vector:         coupon vector 
%
% OUTPUT:
% swaptionPrices:   swaption prices          
%
% USES
% disc_factor_fwd_swaption
%

% HW param a 
a = jData.a;
% HW param sigma  
sigma = jData.sigma;

% Coupon bond put strike 
strike = jData.strike; 

% Yearfrac from t0 to every payment dates with ACT 365
delta_YF365 = jData.delta_YF365;


% Spot discount factor vector
spot_discount = jData.spot_discount;

% Coupon vector
c_vector = jData.c_vector;

% Present the swaption prices vector from 2y10y to 9y10y
swaptionPrices=zeros(8,1);

for start_year = 2:9
    % Compute a yearfrac vector with needed delta
    delta_y_365=delta_YF365(start_year:end);
    
    % Compute a vector with needed coupon 
    c_swaption = c_vector((start_year-2)+1:end);
    
    % Forward discount factor in t0 from starting year to 10y
    disc_y = disc_factor_fwd_swaption(spot_discount,start_year);
    
    % Hull-White sigma
    sigma_HW = @(s,t) sigma.*(1-exp(-a*(t-s)))./a;
    
    % Evaluate x_star such that P(x_star) = K
    
    % Forward discount factors with first date t  
    B = @(x) disc_y.*exp(-x.*sigma_HW(0,delta_y_365(2:end)-delta_y_365(1))/sigma-...
              0.5.*integral(@(u) sigma_HW(u,delta_y_365(2:end)).^2 - sigma_HW(u,delta_y_365(1)).^2,0,delta_y_365(1),'ArrayValued',true)); 
    % Particular underlying in Jamshidian approach
    P = @(x) c_swaption'*B(x);
    % Evaluate the x_star with fzero function 
    x_star = fzero (@(x) P(x) - strike , 0);
    % Compute strikes
    B_star= B(x_star);
    
    % Spot discount (t0,starting year)
    spot_y = spot_discount(start_year);
     
    % Evaluate the call with Gaussian HJM formula
    V = sqrt(1/delta_y_365(1).*integral(@(u) (sigma_HW(u,delta_y_365(2:end)) - sigma_HW(u,delta_y_365(1))).^2,0,delta_y_365(1),'ArrayValued',true));
    d1 = log (disc_y./B_star)./(V.*sqrt(delta_y_365(1))) + 0.5.*V.*sqrt(delta_y_365(1));
    d2 = log (disc_y./B_star)./(V.*sqrt(delta_y_365(1))) - 0.5.*V.*sqrt(delta_y_365(1));
    ZCB = spot_y*(disc_y.*normcdf(d1) - B_star.*normcdf(d2)); 
    
    % Compute the price of coupon bond call
    c_call = c_swaption'*ZCB;
    
    % Put-Call parity to have the swaption price
    swaptionPrices(start_year-1) =c_call- c_swaption'*spot_discount(start_year+1:end) + spot_discount(start_year);
end

end %swaptionPriceJamshidian