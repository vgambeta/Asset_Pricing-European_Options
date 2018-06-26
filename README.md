# Asset_Pricing-European_Options
A comparison between monte-carlo pricing and balck-scholes model pricing of a non-dividend paying European option

The purpose of this task is to compare the results of pricing a non-dividend paying European options based on two different models. The analytical Black-Scholes pricing model and the simulated Monte-Carlo pricing model are used to find the expected price of the options and is implemented in MATLAB.

Black-Scholes Pricing Model:

The Black-Scholes model is widely used and is solved analytically. A direct computation of the prices can be done using this model. This is implemented in MATLAB and is very straight forward. It is considered the based model to use and will be the base line comparison price to the Monte Carlo simulation method.

Monte-Carlo Pricing Model:

Using Geometric Brownian Motion (GBM) with a constant drift (μ) and volatility (𝜎𝜎) with Monte Carlo simulations can be used to price the underlying stock and then computing the payoffs of the options and finish by discounting the payoffs back to get the price of the option. 

The random walk paths are generated using a provided function which uses the above formulation to generate the price each step. This set of paths are then used to calculate the pay-offs which are discounted back using the risk-free rate to determine the price of the option.

Barrier options are also explored using a Monte-Carlo simualtion.

Results:
The Black-Scholes method is much less computationally involved and produces similar results to the simulation. All three methods are within a tight tolerance of each other which and is considered realistic.
