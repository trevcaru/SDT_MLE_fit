#!/usr/bin/env python
# coding: utf-8

# This code is a branch from the original package. Ported to Python by Trevor Caruso. Please contact trevorcaruso@ufl.edu for any questions or comments.
# 
# For more information about metacognitive and Type 2 Signal Detection Theory (SDT), visit the following website at http://www.columbia.edu/~bsm2105/type2sdt/ and refer to the following publications:
# 
# Maniscalco, B., & Lau, H. (2012). "A signal detection theoretic approach for estimating metacognitive sensitivity from confidence ratings." Published in Consciousness and Cognition, 21(1), 422–430. doi:10.1016/j.concog.2011.09.021
# 
# Maniscalco, B., & Lau, H. (2014). "Signal detection theory analysis of type 1 and type 2 data: meta-d’, response-specific meta-d’, and the unequal variance SDT mode." In S. M. Fleming & C. D. Frith (Eds.), The Cognitive Neuroscience of Metacognition (pp.25-66). Published by Springer.
# 
# If you utilize these functions, please cite the aforementioned papers and scripts upon which they are based.
# 
# The same input-format as M & L is expected. As per the original code:

# In[ ]:


# INPUTS
#
# * nR_S1, nR_S2
# these are vectors containing the total number of responses in
# each response category, conditional on presentation of S1 and S2.
#
# e.g. if nR_S1 = [100 50 20 10 5 1], then when stimulus S1 was
# presented, the subject had the following response counts:
# responded S1, rating=3 : 100 times
# responded S1, rating=2 : 50 times
# responded S1, rating=1 : 20 times
# responded S2, rating=1 : 10 times
# responded S2, rating=2 : 5 times
# responded S2, rating=3 : 1 time
#
# The ordering of response / rating counts for S2 should be the same as it
# is for S1. e.g. if nR_S2 = [3 7 8 12 27 89], then when stimulus S2 was
# presented, the subject had the following response counts:
# responded S1, rating=3 : 3 times
# responded S1, rating=2 : 7 times
# responded S1, rating=1 : 8 times
# responded S2, rating=1 : 12 times
# responded S2, rating=2 : 27 times
# responded S2, rating=3 : 89 times


# The function's input should consist of counts for each of these responses separately for each stimulus type.

# In[ ]:


import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

def SDT_MLE_fit(nR_S1, nR_S2):
    nRatings = len(nR_S1) // 2
    nCriteria = 2 * nRatings - 1

    nR_S1 = np.array(nR_S1)
    nR_S2 = np.array(nR_S2)

    ratingHR = np.cumsum(nR_S2[::-1]) / np.sum(nR_S2)
    ratingFAR = np.cumsum(nR_S1[::-1]) / np.sum(nR_S1)
    s = np.polyfit(norm.ppf(ratingFAR[:-1]), norm.ppf(ratingHR[:-1]), 1)[0]

    d1 = np.mean((1 / s) * norm.ppf(ratingHR[:-1]) - norm.ppf(ratingFAR[:-1]))
    c1 = np.linspace(-3, 3, nCriteria)  
    guess = np.concatenate(([d1, s], c1))

    def SDT_logL(parameters):
        # Parameter values must be within realistic bounds to avoid 'nan' likelihood
        d, s = parameters[0], parameters[1]
        if s <= 0 or d <= 0:
            return np.inf  # return infinity to denote invalid log likelihood
        
        S1mu = -parameters[0] / 2
        S2mu = parameters[0] / 2
        S1sd = 1
        S2sd = 1 / parameters[1]
        c1 = parameters[2:]
        
        ci = np.concatenate(([-np.inf], c1))
        cj = np.concatenate((c1, [np.inf]))
        
        pC_S1 = norm.cdf(cj, S1mu, S1sd) - norm.cdf(ci, S1mu, S1sd)
        pC_S2 = norm.cdf(cj, S2mu, S2sd) - norm.cdf(ci, S2mu, S2sd)
        
        # Safeguard against taking log of zero
        pC_S1 = norm.cdf(cj, S1mu, S1sd) - norm.cdf(ci, S1mu, S1sd)
        pC_S2 = norm.cdf(cj, S2mu, S2sd) - norm.cdf(ci, S2mu, S2sd)
        
        logL = np.sum(nR_S1 * np.log(np.maximum(pC_S1, 1e-10))) + np.sum(nR_S2 * np.log(np.maximum(pC_S2, 1e-10)))
        return -logL
        
    print(s)
    # Define the bounds for the parameters
    bounds = [(1e-5, None), # bounds for d
              (1e-5, None)] # bounds for s
    bounds += [(-3, 3) for _ in range(nCriteria)] # bounds for c1

    result = minimize(SDT_logL, guess, method='L-BFGS-B', bounds=bounds, options={'maxiter': 10000})
    

    # Initialize result_dict outside the if-else block
    result_dict = {
        's': None,
        'logL': None,
        'k': nCriteria,
        'n': np.sum(nR_S1) + np.sum(nR_S2),
        'AIC': None,
        'BIC': None,
        'success': result.success,
        'message': result.message
    }

    if not result.success:
        print(f"Optimization failed: {result.message}")
    else:
        optimized_params = result.x
        S1mu = -optimized_params[0] / 2
        S2mu = optimized_params[0] / 2
        S1sd = 1
        S2sd = 1 / optimized_params[1]
        c1 = optimized_params[2:]
        ci = np.concatenate(([-np.inf], c1))
        cj = np.concatenate((c1, [np.inf]))
        pC_S1 = norm.cdf(cj, S1mu, S1sd) - norm.cdf(ci, S1mu, S1sd)
        pC_S2 = norm.cdf(cj, S2mu, S2sd) - norm.cdf(ci, S2mu, S2sd)
        pC_S1 = np.maximum(pC_S1, 1e-100)
        pC_S2 = np.maximum(pC_S2, 1e-100)
        logL = np.sum(nR_S1 * np.log(pC_S1)) + np.sum(nR_S2 * np.log(pC_S2))
        
        result_dict['s'] = optimized_params[1]
        result_dict['logL'] = -logL
        n = result_dict['n']
        AIC = 2 * nCriteria - 2 * logL
        BIC = nCriteria * np.log(n) - 2 * logL
        
        result_dict['AIC'] = AIC
        result_dict['BIC'] = BIC
        print(f"s value: {optimized_params[1]}")
        print(f"logL: {-logL}")
        print(f"AIC: {AIC}")
        print(f"BIC: {BIC}")

    return result_dict


# In[20]:


# Example usage with the provided nR_S1 and nR_S2 (sanity check; result should be ~1)
nR_S1 = [100, 50, 20, 10, 5, 1]
nR_S2 = [1, 5, 10, 20, 50, 100]

result = SDT_MLE_fit(nR_S1, nR_S2)
print(result)


# If using the code to extract the estimated s-value:

# In[21]:


# Extracting the s-value from the result and storing it in the variable s_estim
s_estim = result['s']
print("s_estim:", s_estim)

