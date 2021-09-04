//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Random Number Distributions

This chapter describes functions for generating random variates and computing their probability distributions.
Samples from the distributions described in this chapter can be obtained using any of the random number generators in the library as an underlying source of randomness.

In the simplest cases a non-uniform distribution can be obtained analytically from the uniform distribution of a random number generator by applying an appropriate transformation.
This method uses one call to the random number generator. More complicated distributions are created by the acceptance-rejection method, which compares the desired distribution against a distribution which is similar and known analytically.
This usually requires several samples from the generator.

The library also provides cumulative distribution functions and inverse cumulative distribution functions, sometimes referred to as quantile functions.
The cumulative distribution functions and their inverses are computed separately for the upper and lower tails of the distribution, allowing full accuracy to be retained for small results.

Note that the discrete random variate functions always return a value of type unsigned int, and on most platforms this has a maximum value of 2^32-1 ~=~ 4.29e9. They should only be called with a safe range of parameters (where there is a negligible probability of a variate exceeding this limit) to prevent incorrect results due to overflow.

## Introduction

Continuous random number distributions are defined by a probability density function, p(x), such that the probability of x occurring in the infinitesimal range x to x+dx is p dx.

The cumulative distribution function for the lower tail P(x) is defined by the integral,

P(x) = \int_{-\infty}^{x} dx' p(x')

and gives the probability of a variate taking a value less than x.

The cumulative distribution function for the upper tail Q(x) is defined by the integral,

Q(x) = \int_{x}^{+\infty} dx' p(x')

and gives the probability of a variate taking a value greater than x.

The upper and lower cumulative distribution functions are related by P(x) + Q(x) = 1 and satisfy 0 <= P(x) <= 1, 0 <= Q(x) <= 1.

The inverse cumulative distributions, x=P^{-1}(P) and x=Q^{-1}(Q) give the values of x which correspond to a specific value of P or Q. They can be used to find confidence limits from probability values.

For discrete distributions the probability of sampling the integer value k is given by p(k), where \sum_k p(k) = 1. The cumulative distribution for the lower tail P(k) of a discrete distribution is defined as,

P(k) = \sum_{i <= k} p(i)

where the sum is over the allowed range of the distribution less than or equal to k.


The cumulative distribution for the upper tail of a discrete distribution Q(k) is defined as

Q(k) = \sum_{i > k} p(i)

giving the sum of probabilities for all values greater than k. These two definitions satisfy the identity P(k)+Q(k)=1.

If the range of the distribution is 1 to n inclusive then P(n)=1, Q(n)=0 while P(1) = p(1), Q(1)=1-p(1).
!*/

pub mod bernoulli;
pub mod beta;
pub mod binomial;
pub mod bivariate_gaussian;
pub mod cauchy;
pub mod chi_squared;
pub mod dirichlet;
pub mod exponential;
pub mod exponential_power;
pub mod f_distribution;
pub mod flat;
pub mod gamma;
pub mod gaussian;
pub mod gaussian_tail;
pub mod geometric;
pub mod gumbel;
pub mod hypergeometric;
pub mod landau;
pub mod laplace;
pub mod logarithmic;
pub mod logistic;
pub mod lognormal;
pub mod multinomial;
pub mod negative_binomial;
pub mod pareto;
pub mod pascal;
pub mod poisson;
pub mod rayleigh;
pub mod rayleigh_tail;
pub mod t_distribution;
pub mod weibull;
