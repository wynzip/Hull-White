function spot_discount = disc_factor_spot(dates, discounts, new_dates)
% Compute spot discount factors for new dates thanks to given informations
%
%INPUT
% dates:            dates vector
% discounts:        discount factors vector corresponding to given dates 
% new_dates:        dates for which we want the discount factors
%
%OUTPUT
% spot_discount:    needed discount factors     
% 

    % date convetion
    C_ACT_365 = 3;
    
    % transforming type for yearfrac
    dates = datenum(dates);
    new_dates = datenum(new_dates);
    
    % obtaining zero rates from bootstrapped data
    year_frac = yearfrac(dates(1), dates(2:end), C_ACT_365);  %delta( t_0, t_market)
    zero_rates= -log(discounts(2:end))./year_frac; %zero rates from market
    
    % interpolating on zero rates
    new_zero_rates = interp1(dates(2:end), zero_rates, new_dates, 'linear'); % zero rates for payment dates
    
    % obtaining discount factors (spot) realtive to payment dates
    year_frac_spot = yearfrac(dates(1), new_dates, C_ACT_365);  %delta( t_0, t_payments)
    
    %discount factors corresponding to dates of payments
    spot_discount = exp(-new_zero_rates .* year_frac_spot); 

end %disc_factor_spot 