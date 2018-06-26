# Asset_Pricing-European_Options
A comparison between monte-carlo pricing and balck-scholes model pricing of a non-dividend paying European option

The purpose of this task is to compare the results of pricing a non-dividend paying European options based on two different models. The analytical Black-Scholes pricing model and the simulated Monte-Carlo pricing model are used to find the expected price of the options and is implemented in MATLAB.

Black-Scholes Pricing Model:

The Black-Scholes model is widely used and is solved analytically. A direct computation of the prices can be done using this model. This is implemented in MATLAB and is very straight forward. It is considered the based model to use and will be the base line comparison price to the Monte Carlo simulation method.

Monte-Carlo Pricing Model:

Using Geometric Brownian Motion (GBM) with a constant drift (μ) and volatility (𝜎𝜎) with Monte Carlo simulations can be used to price the underlying stock and then computing the payoffs of the options and finish by discounting the payoffs back to get the price of the option. 

The random walk paths are generated using a provided function which uses the formulation to generate the price at each step. This set of paths are then used to calculate the pay-offs which are discounted back using the risk-free rate to determine the price of the option.

Barrier options are also explored using a Monte-Carlo simualtion.

Results:
The Black-Scholes method is much less computationally involved and produces similar results to the simulation. All three methods are within a tight tolerance of each other which and is considered realistic.


%FUNCTIONS
function [BS_CALL_PRICE, BS_PUT_PRICE] = BS_european_price(S0, K, T, r, sigma)
%Compute Parameters
D1 = ( 1/(sigma*sqrt(T)) ) * ( log( S0./K ) + (r + (sigma^2)/2) * T);
D2 = D1 - sigma .* sqrt(T);
%Compute Call & Put Price
BS_CALL_PRICE = S0.*normcdf(D1,0,1) - K.*exp(-r.*T).*normcdf(D2,0,1);
BS_PUT_PRICE = K.*exp(-r.*T).*normcdf(-D2,0,1) - S0.*normcdf(-D1,0,1);
%Built in BLS Function (Confirms)
%[Call, Put] = blsprice(100, 105, 0.05, 1, 0.2)
end
function [MC_CALL_PRICE, MC_PUT_PRICE] = MC_european_price(S0, K, T, r, mu, sigma, numPaths)
%Initialize Steps
numSteps = 1;
%Generate GRW Paths using Function
MC_PATH = GRWPaths(S0, mu, sigma, T, numSteps, numPaths);
%Calculate Pay-Off of Option for Each Price Path
CALL_PAY = max(MC_PATH(end,:) - K, 0);
PUT_PAY = max(K - MC_PATH(end,:), 0);
%Dicsount Pay-Off Back to Time 0 to get Price of Option
MC_CALL_PRICE = mean(CALL_PAY) * exp(-r.*T);
MC_PUT_PRICE = mean(PUT_PAY)* exp(-r.*T);
end
function [MC_CALL_PRICE, MC_PUT_PRICE] = MC_european_price_multi(S0, K, T, r, mu, sigma, numSteps, numPaths)
%Generate GRW Paths using Function
MC_PATH = GRWPaths(S0, mu, sigma, T, numSteps, numPaths);
%Calculate Pay-Off of Option for Each Price Path
CALL_PAY = max(MC_PATH(end,:) - K, 0);
PUT_PAY = max(K - MC_PATH(end,:), 0);
%Dicsount Pay-Off Back to Time 0 to get Price of Option
MC_CALL_PRICE = mean(CALL_PAY) * exp(-r.*T);
MC_PUT_PRICE = mean(PUT_PAY)* exp(-r.*T);
end
function [MC_CALL_PRICE, MC_PUT_PRICE] = MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths)
%Adjust Volatility
%sigma = sigma - 0.1;
%Generate GRW Paths using Function
MC_PATH = GRWPaths(S0, mu, sigma, T, numSteps, numPaths);
%Check Each Path for Boundary Cross (Knock-In)
[rows, columns] = ind2sub(size(MC_PATH), find(MC_PATH >= Sb));
%Find Columns that Need Set to 0 (Difference between)
idx = setdiff([1:1:numPaths], unique(columns));
%Keep Paths that Cross Boundary
MC_PATH(:,idx) = 0;
%Remove Zeros from Simulated Price List
PATHS = nonzeros(MC_PATH(end,:));
%Calculate Pay-Off of Option for Each Price Path
CALL_PAY = max(PATHS - K, 0);
PUT_PAY = max(K - PATHS, 0);
%Add Back Zeros to Correct Array Length for Average Calculation
CALL_PAY = [CALL_PAY; zeros((numPaths - length(PATHS)),1)];
PUT_PAY = [PUT_PAY; zeros((numPaths - length(PATHS)),1)];
%Dicsount Pay-Off Back to Time 0 to get Price of Option
MC_CALL_PRICE = mean(CALL_PAY) * exp(-r.*T);
MC_PUT_PRICE = mean(PUT_PAY) * exp(-r.*T);
End
Appendix B – Monte Carlo Algorithm to Match Black-Scholes
while abs(callMC_N_AVG - BS_CALL) > 0.01
%for i = 1:numScenarios
[callMC_European_Price_multi_step, putMC_European_Price_multi_step, multi_paths, call_pay, put_pay] = MC_european_price_multi(S0, K, T, r, mu, sigma, numSteps, numPaths);
callMC_N(cnt) = callMC_European_Price_multi_step;
putMC_N(cnt) = putMC_European_Price_multi_step;
cnt = cnt +1;
%Find Mean Call & Put Price of N Scenarios
callMC_N_AVG = mean(callMC_N);
putMC_N_AVG = mean(putMC_N);
end
