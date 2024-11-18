function discount_fwrd = disc_factor_fwd_swaption(factor_spot,start_year)
% Evaluate forward discount factor from fixed year to 10y
%
% INPUT
% factor_spot      spot discount factors  
% start_year       start year 
%
% OUTPUT
% discount_fwrd    needed forward discount factors   
%

discount_fwrd = factor_spot(start_year+1:end)/factor_spot(start_year);

end %disc_factor_fwd_swaption