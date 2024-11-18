function upfront = computeNPV_3y(mktData, swapData, p)
% This function computes the upfront from the npv of the swap described 
% in Annex 2, for contract with maturity 3 years.
%
% INPUT:
% mktData:          struct with needed discounts and yearfracs
% swapData:         struct with needed values of the swap contract
% p:                vector with probability of all possible scenarios
%
% OUTPUT:
% upfront:          value of upfront of the 3y swap

    % Prepare the variables: discounts and yearfracs
    B = mktData.payDiscounts;
    yf_360 = mktData.payYF_360;
    yearly_yf = mktData.couponYF_30_360;

    % Prepare the variables: spol, coupon1=6%, coupon2=2%
    spol = swapData.spolA;
    c6 = swapData.coupon1;
    c2 = swapData.coupon2;

    % Compute inter time yeafracs
    yf_1y = [yearly_yf(1); yearly_yf(2:end) - yearly_yf(1:end-1)];
    yf_3m = [yf_360(1); yf_360(2:end) - yf_360(1:end-1)];

    % Convention: leg B (the upfront one) has a positive value
    % Compute npv in the up_up case
    npv_up_up = @(x) x + c2*yf_1y(3)*B(12) - 1 + B(12) - spol*yf_3m'*B;

    % Compute npv in the up_down case
    npv_up_down = @(x) x + c6*yf_1y(2)*B(8) - 1 + B(8) - spol*yf_3m(1:8)'*B(1:8);

    % Compute npv in the down_up case
    npv_down_up = @(x) x + c6*yf_1y(1)*B(4) - 1 + B(8) - spol*yf_3m(1:8)'*B(1:8);

    % Compute npv in the down_down case
    npv_down_down = @(x) x + c6*yf_1y(1)*B(4) + c6*yf_1y(1)*B(8) - 1 + B(8) ...
        - spol*yf_3m(1:8)'*B(1:8); 

    % Compute npv with all the probabilities
    npv_total = @(x) p(4)*npv_up_up(x) + p(3)*npv_up_down(x) + p(2)*npv_down_up(x) ...
        + p(1)*npv_down_down(x);
    
    % Find the upfront X
    upfront = fzero(npv_total,0.033);

end