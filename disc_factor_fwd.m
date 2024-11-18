function fwd_discounts = disc_factor_fwd(dates, discounts, year_frac)
% Compute fwd discount factors for new dates thanks to given informations
%
%INPUT
% dates:            dates vector
% discounts:        discount factors vector corresponding to given dates 
% new_dates:        dates for which we want the discount factors
%
%OUTPUT
% spot_discount:    needed discount factors     
% 

    C_ACT_365 = 3;

    % obtaining zero rates from bootstrapped data
    yf = yearfrac(dates(1), dates(2:end), C_ACT_365);  %delta( t_0, t_market)
    zero_rates= -log(discounts(2:end))./yf; %zero rates from market
    
    % interpolating on zero rates
    new_zero_rates = interp1(yf, zero_rates, year_frac, 'linear'); % zero rates for payment dates
    
    %discount factors corresponding to dates of payments
    spot_discount = exp(-new_zero_rates .* year_frac); 

    % compute fwd discounts
    fwd_discounts = spot_discount(2:end)./spot_discount(1:end-1);
    fwd_discounts = [spot_discount(1); fwd_discounts];

end %disc_factor_spot 