clc 
clear
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           EXERCISE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('...Exercise 1...')
%%
% load data of stock
load('cSelect20230131_B.mat');

% load bootstrap data
load('dates.mat');
load('discounts.mat');
bootstrap.dates = dates;
bootstrap.discounts = discounts;
ACT_360 = 2;
ACT_365 = 3;
C_30_360 = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Present data

% Time data
TTM = 3; % Time to maturity in years
payments_step = 3; % Payments step in months
dateSet = dates(1); % settlement date
schedule_3m = findDates(dateSet,payments_step,TTM,0); % 3-monthly payments dates, including current date
% Reset dates (2 days before the yearly payment dates)
dateReset = datetime(dates(1),'ConvertFrom','datenum') + calyears(1:TTM) -2; 

% Set up mkt Data struct
% spot discount factor vector in payment dates
mktData.payDiscounts = disc_factor_spot(dates, discounts, schedule_3m(2:end));
% discounts on reset dates
mktData.resetDiscounts = disc_factor_spot(dates, discounts, dateReset);
% year fractions
mktData.payYF_360 = yearfrac(dates(1), schedule_3m(2:end), ACT_360);
mktData.payYF_30_360 = yearfrac(dates(1), schedule_3m(2:end), C_30_360);
mktData.couponYF_30_360 = yearfrac(dateSet, schedule_3m([5,9,13]),C_30_360);
mktData.resetYF_365 = yearfrac(dates(1),dateReset,ACT_365)';

% Set up Swap Data struct, info about our structured product
swapData.spolA = 0.013;
swapData.paymentDatesA = schedule_3m(2:end);
%swapData.paymentDatesB = findDates(dateSet,12,TTM,0);
swapData.settlementDate = dates(1);
swapData.Notional = 100e6;
swapData.coupon1 = 0.06;
swapData.coupon2 = 0.02;
swapData.strike = 3200; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calibrate NIG and VG model parameters

% Obtaining the needed data from the given struct
d = cSelect.dividend; % dividend
K_mkt = cSelect.reference; % strike
S0 = K_mkt; % spot price (ATM Spot condition)

% Finding the istantaneous interest rate r for the reset date in 1 year (= zero rate)
TT_reset1 = yearfrac(dates(1),dateReset(1), 3); % time to first reset
zRates = zeroRates(dates,discounts); % find whole zeroRates curve from bootstrap
r = interp1(yearfrac(dates(1),dates(2:end),3),zRates,TT_reset1)/100; % linear interpolation

% Setting up the inputs of the FFT pricing method we will use in the
% lsqnonlin least squares optimization method
inputStruct.deltaT = 1; % denominator of variance of G from text
inputStruct.B = exp(-r*TT_reset1);  % Discount factor B(t0,t) where t0=dateSet and t=dateMat 
inputStruct.F0 = S0*exp(-d*TT_reset1)/inputStruct.B;  % Forward at t0

% Log-moneyness grid where we evaluate the Call prices
x = log(inputStruct.F0./cSelect.strikes)'; 

% Computation of market prices for Call following Black formula
mktCall = blkprice(inputStruct.F0, cSelect.strikes, r, TT_reset1, cSelect.surface)';

% Compute model parameters NIG VG
[paramsNIG,paramsVG] = calibrate_NIG_VG_params(inputStruct,mktCall,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute probability of St < K at 1y with Lewis fromula

% Prepare parameters
deltaT = mktData.resetYF_365(1);
TTM = 1; % fixed by hand
K_swap = 3200;
kappa = log(S0/K_swap) + (r-d)*TTM;

% Compute probability for NIG model with Lewis
p_down_1y_NIG_Lewis = pricingLewisProb(kappa,deltaT,paramsNIG,1)

% Compute probability for NIG model with Lewis
p_down_1y_VG_Lewis = pricingLewisProb(kappa,deltaT,paramsVG,2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evaluate probability of St < K at 1y with MC simulation of underlying

% Prepare parameters
F0 = inputStruct.F0;
nSim = 1e6;

% Evaluate probability of St < K at 1y with MC simulation of NIG underlying
[p_down_NIG_1y,IC_prob_down_NIG_1y] = simul_MC_NIG_2y(paramsNIG,F0,K_swap,nSim,mktData.resetYF_365(1))

% Evaluate probability of all 4 scenarios at 3y with MC sim of NIG underlying
[p_down_down,p_down_up,p_up_down,p_up_up] = simul_MC_NIG_3y(paramsNIG,F0,K_swap,nSim,mktData.resetYF_365)

% Evaluate probability of St < K at 1y with MC simulation of VG underlying
[p_down_VG_1y,IC_prob_down_VG_1y] = simul_MC_VG(paramsVG,F0,K_swap,nSim,mktData.resetYF_365(1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% a. Compute upfront with NIG model (time to maturity is 2 years)

% Value of upfront (all computed inside the compute upfront function) 
upfront_MC_NIG_2y = computeNPV_2y(mktData, swapData, p_down_NIG_1y)

% Compute confidence interval for this upfront
lb = computeNPV_2y(mktData, swapData, IC_prob_down_NIG_1y(2));
ub = computeNPV_2y(mktData, swapData, IC_prob_down_NIG_1y(1));
IC_upfront_MC_NIG_2y = [lb; ub]

%consider the principal amount
total_upfront_MC_NIG_2y = upfront_MC_NIG_2y* swapData.Notional;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% b. Compute upfront with VG model (time to maturity is 2 years)

% Value of upfront (all computed inside the compute upfront function) 
upfront_MC_VG_2y = computeNPV_2y(mktData, swapData, p_down_VG_1y)

% Compute confidence interval for this upfront
lb = computeNPV_2y(mktData, swapData, IC_prob_down_VG_1y(2));
ub = computeNPV_2y(mktData, swapData, IC_prob_down_VG_1y(1));
IC_upfront_MC_VG_2y = [lb; ub]

%consider the principal amount
total_upfront_MC_VG_2y = upfront_MC_VG_2y*swapData.Notional;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% d. Compute upfront with NIG model (time to maturity is 3 years)

% Value of upfront (all computed inside the compute upfront function) 
upfront_MC_3y = computeNPV_3y(mktData, swapData, [p_down_down,p_down_up,p_up_down,p_up_up])

% consider the principal amount
total_upfront_MC_3y = upfront_MC_3y*swapData.Notional;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% e. compute upfront using black formula (time to maturity is 2 years)

% Compute upfront of the swap under the assumption of Black model
upfront_blk = compute_upfront_blk(mktData,swapData,cSelect,F0,dates,dateReset,zRates)

% consider the principal amount
total_upfront_blk = upfront_blk*swapData.Notional;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           EXERCISE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('...Exercise 2...')
%% a. Compute Bermudan Swaption price through a trinomial HW tree
% Let's clean the workspace 
clear
load('dates.mat')
load('discounts.mat')

% Number of the time step in the tree from 0 to 10 years
M = 3650;
% Parameters of Hull and White model
a = 0.11;
sigma = 0.008;

% Initialize needed parameters of the tree
[x_tree,stoc_disc_mat,yearly_fwd_disc,dt,tree_yearfracs,yf_365_yearly] = initialize_Tree(M,dates,discounts,a,sigma);

% Mu_hat for the OU process
mu_hat = 1 - exp(-a*dt);
% Recompute bounds of the tree
l_max = (1-sqrt(2/3))/mu_hat;
l_max = ceil(l_max);
years = 10;

% Compute bermudan swaption price thanks to the tree
price_bermudan = bermudan_tree(x_tree,stoc_disc_mat,yearly_fwd_disc,dt,tree_yearfracs,yf_365_yearly,mu_hat,l_max,years,M,dates,a,sigma)

% Plot for the convergence of the Tree
%convergence_tree(dates,discounts)
%this line is commented as it takes a long time  

%% check on the european price for a fwd starting swap (9-10)
price_european = european_tree(x_tree,stoc_disc_mat,yearly_fwd_disc,dt,tree_yearfracs,mu_hat,l_max,M,dates,a,sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% c. Jamshidian formula
bootstrap.dates = dates;
bootstrap.discounts = discounts;
ACT_360 = 2;
ACT_365 = 3;
ACT_30_360 = 6;

% Create Struct for Jamshidian approach
% HW param a 
jData.a = 0.11;
% HW param sigma  
jData.sigma = 0.008;  
% Coupon bond put strike 
jData.strike = 1; 
% Swaption strike
jData.c = 0.05; 
% Payment dates annually untill 10th year   
jData.payment_date = findDates(dates(1),12,10,0); 

% Compute different yearfracs

% yearfrac from t0 to every payment dates with ACT 365
jData.delta_YF365 = yearfrac(jData.payment_date(1),jData.payment_date(2:end),ACT_365);

% yearfrac between two consecutive payment dates with ACT 30/360
interdelta_YF30360 = yearfrac(jData.payment_date(1:end-1),jData.payment_date(2:end),ACT_30_360); 
jData.interdelta_YF30360 = interdelta_YF30360(3:end); 

% yearfrac between two consecutive payment dates with ACT 365
interdelta_YF365 = yearfrac(jData.payment_date(1:end-1),jData.payment_date(2:end),ACT_365); 
jData.interdelta_YF365 = interdelta_YF365(3:end);  

% Spot discount factor vector from t0 and payment dates
jData.spot_discount = disc_factor_spot(dates, discounts, jData.payment_date(2:end));

% Coupon vector  
c_vector = jData.c*ones(length(jData.interdelta_YF30360),1);
c_vector = c_vector.* jData.interdelta_YF30360;
c_vector(end) = c_vector(end) + 1;
jData.c_vector = c_vector;

% Swaption prices with Jamshidian approach
prices_swaption = swaptionPriceJamshidian(jData);


% Evaluation of bounds
upper_bound = sum(prices_swaption)
lower_bound = max(prices_swaption)



