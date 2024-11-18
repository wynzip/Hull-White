function upfront = computeNPV_2y(mktData, swapData, p)
% This function computes the upfront from the npv of the swap described 
% in Annex 2, for contract with maturity 2 years.
%
% INPUT:
% mktData:          struct with needed discounts and yearfracs
% swapData:         struct with needed values of the swap contract
% p:                vector with probability of all possible scenarios
%
% OUTPUT:
% upfront:          value of upfront of the 2y swap

    % NPV of payment leg made by party A
    legA = compute_leg_A_2y(mktData, swapData, p);
    % NPV of payment leg made by party B
    legB = compute_leg_B_2y(mktData, swapData, p);
    
    % Value of upfront 
    upfront = legA - legB;

end