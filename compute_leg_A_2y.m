function legA = compute_leg_A_2y(mktData, swapData, probs)
% This function computes the NPV of the payment leg made by party A in the
% swap described in Annex 2, for contract with maturity 2 years.
%
% INPUT:
% mktData:          struct with needed discounts and yearfracs
% swapData:         struct with needed values of the swap contract
% probs:            vector with probability of stock<strike for considered
%                   years, according to all possible scenarios
%
% OUTPUT:
% legA:             value of leg A of the contract


    %unpack needed data
    spolA = swapData.spolA;
    interFracs = [mktData.payYF_360(1); mktData.payYF_360(2:end) - mktData.payYF_360(1:end-1)];
    payDisc = mktData.payDiscounts; % a value every 3 months

    % number of flows in a year
    nFlows = 4;
    idx_1y = nFlows;
    idx_2y = 2*nFlows;
            
    % Compute leg A
    legA = probs(1)*(spolA*interFracs(1:idx_1y)'*payDisc(1:idx_1y) + 1 - payDisc(idx_1y)) +...
         (1-probs(1))*(spolA*interFracs(1:idx_2y)'*payDisc(1:idx_2y) + 1 - payDisc(idx_2y));
  

end % compute_leg_A