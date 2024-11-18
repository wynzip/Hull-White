function upfront_blk = compute_upfront_blk(mktData,swapData,cSelect,F0,dates,dateReset,zRates)
% This function computes the upfront with Black model of the swap described 
% in Annex 2, for contract with maturity 2 years.
%
% INPUT:
% mktData:          struct with needed discounts and yearfracs
% swapData:         struct with needed values of the swap contract
% cSelect:          structure with some market data for the Black model      
% F0:               fwd value in 0
% dates:            dates of the bootstrap
% dateReset:        coupon reset dates (of the swap)
% zRates:           zero rates from the bootstrap
%
% OUTPUT:
% upfront:          value of upfront of the 2y swap


    % epsilon small step
    eps = 10;
    % compute mkt volatilities that correspond to the strikes considered for
    % the digital
    strikes_dig_blk = [swapData.strike, swapData.strike+eps];
    vols_dig_blk = interp1(cSelect.strikes, cSelect.surface, strikes_dig_blk, "spline");
    
    % time to first reset dates, and corresponding zero rate
    TT_reset1 = yearfrac(dates(1),dateReset(1), 3); % time to first reset
    r = interp1(yearfrac(dates(1),dates(2:end),3),zRates,TT_reset1)/100; % linear interpolation
    
    % Find prices of the two Black calls
    calls1_blk = blkprice(F0 , strikes_dig_blk(1), r, TT_reset1,  vols_dig_blk(1));
    calls2_blk = blkprice(F0 , strikes_dig_blk(2), r, TT_reset1,  vols_dig_blk(2));
    
    % Discretized value of the digital option --> P(S_T > K)
    price_digital_blk = (calls1_blk-calls2_blk)/eps;
    
    % Compute the two legs of the contract
    legA_blk = compute_leg_A_2y(mktData, swapData, 1-price_digital_blk);
    legB_blk = compute_leg_B_2y(mktData, swapData, 1-price_digital_blk);
    
    % Compute upfront with Black model
    upfront_blk = legA_blk - legB_blk

end