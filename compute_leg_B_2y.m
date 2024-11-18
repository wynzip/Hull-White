function legB = compute_leg_B_2y(mktData, swapData, probs)
% This function computes the NPV of the payment leg made by party B in the
% swap described in Annex 2, for contract with maturity 2 years.
%
% INPUT:
% mktData:          struct with needed discounts and yearfracs
% swapData:         struct with needed values of the swap contract
% probs:            vector with probability of stock<strike for considered
%                   years, according to all possible scenarios
%
% OUTPUT:
% legB:             value of leg B of the contract

    % number of flows in a year
    nFlows = 4;
    idx_1y = nFlows;
    idx_2y = 2*nFlows;
    
    %unpack needed data
    payDisc = mktData.payDiscounts;
    coupon1 = swapData.coupon1;
    coupon2 = swapData.coupon2;
    yf1 = mktData.payYF_30_360(idx_1y); % from 0 to first reset date
    yf2 = mktData.payYF_30_360(idx_2y) - yf1; % from first to second reset date
        
    % Compute leg B
    legB = probs(1)*coupon1*yf1*payDisc(idx_1y) + (1-probs(1))*coupon2*yf2*payDisc(idx_2y);

end % compute_leg_B