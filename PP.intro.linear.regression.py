# LINEAR REGRESSION

import numpy as np
import matplotlib.pyplot as plt

# Initialize RNG
np.random.seed(123)

# True parameter values
alpha, sigma = 1, 1
beta = [1, 2.5]

# Size of the dataset
size = 100

# Predictor variable
X1 = np.linspace(0, 1, size)
X2 = np.linspace(0, .2, size)

# Simulate outcome variable
Y = alpha + beta[0] * X1 + beta[1] * X2 + np.random.randn(size) * sigma

# Model implementation in PyMC3
from pymc3 import Model, Normal, HalfNormal

# Model object: container for the model random variables
basic_model = Model()

# Context manager: PyMC3 objects are added to the model
with basic_model:

    # Stochastic random variables for the unknown model parameters: with Normal prior distributions for the regression
    # coefficients and half-normal prior distribution for the sd of the observations
    alpha = Normal('alpha', mu = 0, sd = 10)
    beta = Normal('beta', mu = 0, sd = 10, shape = 2)
    sigma = HalfNormal('sigma', sd = 1)

    # Deterministic random variable: expected value of outcome
    mu = alpha + beta[0] * X1 + beta[1] * X2

    # Observed stochastic variable: data LH of the model or sampling distribution of observations
    Y_obs = Normal('Y_obs', mu = mu, sd = sigma, observed = Y)


# Obtaining posterior estimates for the unknown variables in the model: MAP and MCMC

# MAP: mode of the posterior distribution (one point):
# Use of BFGS optimization algorithm
from pymc3 import find_MAP

map_estimate = find_MAP(model = basic_model)

print(map_estimate)

# Use of Powell optimization algorithm
from scipy import optimize

map_estimate = find_MAP(model = basic_model, fmin = optimize.fmin_powell)

print(map_estimate)


# MCMC: gradient-based sampling methods (NUTS)
from pymc3 import NUTS, sample

with basic_model:

    # Obtain starting values via MAP
    start = find_MAP(fmin = optimize.fmin_powell)

    # Instantiate sampler
    step = NUTS(scaling = start)

    # Draw 2000 posterior samples
    trace = sample(2000, step, start = start)

print(trace['alpha'][-5:])

# Posterior analysis of the sampling output

# Simple posterior plot
from pymc3 import traceplot

traceplot(trace)






























