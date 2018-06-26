clc;
clear all;
format long

% Pricing a European option using Black-Scholes formula and Monte Carlo simulations
% Pricing a Barrier option using Monte Carlo simulations

S0 = 100;     % spot price of the underlying stock today
K = 105;      % strike at expiry
mu = 0.05;    % expected return
sigma = 0.2;  % volatility
r = 0.05;     % risk-free rate
T = 1.0;      % years to expiry
Sb = 110;     % barrier

% Define variable numSteps to be the number of steps for multi-step MC
% numPaths - number of sample paths used in simulations

% Implement your Black-Scholes pricing formula (Non-Dividend Paying)
[call_BS_European_Price, putBS_European_Price] = BS_european_price(S0, K, T, r, sigma);

% Implement your one-step Monte Carlo pricing procedure for European option
numPaths = 100000;
[callMC_European_Price_1_step, putMC_European_Price_1_step] = MC_european_price(S0, K, T, r, mu, sigma, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for European option
numSteps = 1; 
numPaths = 100000;
[callMC_European_Price_multi_step, putMC_European_Price_multi_step, multi_paths, call_pay, put_pay] = MC_european_price_multi(S0, K, T, r, mu, sigma, numSteps, numPaths);

%Plot paths
figure;
[frequencyCountsMC1, binLocationsMC1] = hist(call_pay, 40);
bar(binLocationsMC1, frequencyCountsMC1,'FaceAlpha', 0.8, 'FaceColor', 'c'); hold on;
[frequencyCountsMC2, binLocationsMC2] = hist(put_pay, 40);
bar(binLocationsMC2, frequencyCountsMC2,'FaceAlpha', 0.8, 'FaceColor', 'y'); hold on;
line([callMC_European_Price_multi_step, callMC_European_Price_multi_step], [0 1000], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
line([putMC_European_Price_multi_step, putMC_European_Price_multi_step], [0 1000], 'Color', 'b', 'LineWidth', 1, 'LineStyle', '--');
hold on;
box off; axis tight; legend boxoff;
title('Simulated Pay-Offs: 365 Steps, 10,000 Scenarios', 'FontWeight', 'bold');
ylabel('Frequency')
xlabel('Pay-Off ($)')
legend('Call Pay-Off Distribution','Put Pay-Off Distribution','Call Price','Put Price')
text(3,1950,'BREAK: 6000 Max','FontSize', 8)
xlim([0 52])
ylim([0 2000])
hold off;

% Implement your multi-step N-Scenarios Monte Carlo pricing procedure for European option
numScenarios = 1000;
numPaths = 1000;
numSteps = 2;
BS_CALL = call_BS_European_Price;
BS_PUT = putBS_European_Price;
callMC_N = zeros(numScenarios,1);
putMC_N = zeros(numScenarios,1);
cnt = 1;
callMC_N_AVG = 100;
putMC_N_AVG = 100;

tic
while ( abs(callMC_N_AVG - BS_CALL) > 0.01 || abs(putMC_N_AVG - BS_PUT) > 0.01 ) 
%for i = 1:numScenarios
    [callMC_European_Price_multi_step, putMC_European_Price_multi_step, multi_paths, call_pay, put_pay] = MC_european_price_multi(S0, K, T, r, mu, sigma, numSteps, numPaths);
    callMC_N(cnt) = callMC_European_Price_multi_step;
    putMC_N(cnt) = putMC_European_Price_multi_step;
    cnt = cnt +1;
    
    %Find Mean Call & Put Price of N Scenarios
    callMC_N_AVG = mean(callMC_N);
    putMC_N_AVG = mean(putMC_N);
end
abs(BS_CALL - callMC_N_AVG)
abs(BS_PUT - putMC_N_AVG)
cnt
toc

% Implement your one-step Monte Carlo pricing procedure for Barrier option
numSteps = 1;
numPaths = 10000;
[callMC_Barrier_Knockin_Price_1_step, putMC_Barrier_Knockin_Price_1_step] = 
MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths);

% Implement your multi-step Monte Carlo pricing procedure for Barrier option
numSteps = 364;
numPaths = 10000;
[callMC_Barrier_Knockin_Price_multi_step, putMC_Barrier_Knockin_Price_multi_step] = ...
MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths);


disp(['Black-Scholes price of an European call option is ',num2str(call_BS_European_Price)])
disp(['Black-Scholes price of an European put option is ',num2str(putBS_European_Price)])
disp(' ');
disp(['One-step MC price of an European call option is ',num2str(callMC_European_Price_1_step)])
disp(['One-step MC price of an European put option is ',num2str(putMC_European_Price_1_step)])
disp(' ');
disp(['Multi-step MC price of an European call option is ',num2str(callMC_European_Price_multi_step)])
disp(['Multi-step MC price of an European put option is ',num2str(putMC_European_Price_multi_step)])
disp(' ')
%disp(['Multi-step N Mean MC price of an European call option is ',num2str(callMC_N_AVG)])
%disp(['Multi-step N Mean MC price of an European put option is ',num2str(putMC_N_AVG)])
%disp(' ')
disp(['One-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_1_step)])
disp(['One-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_1_step)])
disp(' ')
disp(['Multi-step MC price of an Barrier call option is ',num2str(callMC_Barrier_Knockin_Price_multi_step)])
disp(['Multi-step MC price of an Barrier put option is ',num2str(putMC_Barrier_Knockin_Price_multi_step)])


%FUNCTIONS
function [BS_CALL_PRICE, BS_PUT_PRICE] = BS_european_price(S0, K, T, r, sigma)
   
    %Compute Parameters
    D1 = ( 1/(sigma*sqrt(T)) ) * ( log( S0./K ) + (r + (sigma^2)/2) * T);
    D2 = D1 - sigma .* sqrt(T);

    %Compute Call & Put Price
    BS_CALL_PRICE =  S0.*normcdf(D1,0,1) - K.*exp(-r.*T).*normcdf(D2,0,1);
    BS_PUT_PRICE  =  K.*exp(-r.*T).*normcdf(-D2,0,1) - S0.*normcdf(-D1,0,1);
   
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
    PUT_PAY  = max(K - MC_PATH(end,:), 0);

    %Dicsount Pay-Off Back to Time 0 to get Price of Option
    MC_CALL_PRICE = mean(CALL_PAY) * exp(-r.*T);
    MC_PUT_PRICE  = mean(PUT_PAY)* exp(-r.*T);
    
end

function [MC_CALL_PRICE, MC_PUT_PRICE, MC_PATH, CALL_PAY, PUT_PAY] = MC_european_price_multi(S0, K, T, r, mu, sigma, numSteps, numPaths)

    %Generate GRW Paths using Function
    MC_PATH = GRWPaths(S0, mu, sigma, T, numSteps, numPaths);

    %Calculate Pay-Off of Option for Each Price Path
    CALL_PAY = max(MC_PATH(end,:) - K, 0);
    PUT_PAY  = max(K - MC_PATH(end,:), 0);

    %Dicsount Pay-Off Back to Time 0 to get Price of Option
    MC_CALL_PRICE = mean(CALL_PAY) * exp(-r.*T);
    MC_PUT_PRICE  = mean(PUT_PAY)* exp(-r.*T);

end

function [MC_CALL_PRICE, MC_PUT_PRICE] = MC_barrier_knockin_price(S0, Sb, K, T, r, mu, sigma, numSteps, numPaths)

    %Adjust Volatility
    %sigma = sigma - 0.1;
    
    %Generate GRW Paths using Function
    MC_PATH = GRWPaths(S0, mu, sigma, T, numSteps, numPaths);
    
    %Check Each Path for Boundary Cross (Knock-In)
    [rows, columns] = ind2sub(size(MC_PATH), find(MC_PATH >= Sb));

    %Find Columns that Need Set to 0 (Difference between 1-numPAths array and Columns > Sb)
    idx = setdiff([1:1:numPaths], unique(columns));
    
    %Keep Paths that Cross Boundary
    MC_PATH(:,idx) = 0;

    %Should be Call 8.01 and Put 2.28 (Up and In)
    
    %Different pay off formula for up and in put and up and in call?
    %http://www.derivativepricing.com/blogpage.asp?id=15
    %http://quantcalc.net/Barrier.html

    %Remove Zeros from Simulated Price List
    PATHS = nonzeros(MC_PATH(end,:));
     
    %Calculate Pay-Off of Option for Each Price Path
    CALL_PAY = max(PATHS - K, 0);
    PUT_PAY  = max(K - PATHS, 0);
     
    %Add Back Zeros to Correct Array Length for Average Calculation
    CALL_PAY = [CALL_PAY; zeros((numPaths - length(PATHS)),1)];
    PUT_PAY  = [PUT_PAY; zeros((numPaths - length(PATHS)),1)];
    
    %Dicsount Pay-Off Back to Time 0 to get Price of Option
    MC_CALL_PRICE = mean(CALL_PAY) * exp(-r.*T);
    MC_PUT_PRICE  = mean(PUT_PAY) * exp(-r.*T);

end