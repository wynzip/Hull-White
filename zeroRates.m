function zRates = zeroRates(dates, discounts)
%zero rates computed from the discount factor curve
%
%INPUT
% dates:       dates form the market
% discounts:   discount factors 
%

zRates= - log(discounts(2:end))./ yearfrac(dates(1), dates(2:end), 3)*100;


end %zeroRates