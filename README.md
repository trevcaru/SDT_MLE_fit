# SDT_fit_Py

This code is a branch from the original package. Ported to Python by Trevor Caruso. Please contact trevorcaruso@ufl.edu for any questions or comments.

For more information about metacognitive and Type 2 Signal Detection Theory (SDT), visit the following website at http://www.columbia.edu/~bsm2105/type2sdt/ and refer to the following publications:

Maniscalco, B., & Lau, H. (2012). "A signal detection theoretic approach for estimating metacognitive sensitivity from confidence ratings." Published in Consciousness and Cognition, 21(1), 422–430. doi:10.1016/j.concog.2011.09.021

Maniscalco, B., & Lau, H. (2014). "Signal detection theory analysis of type 1 and type 2 data: meta-d’, response-specific meta-d’, and the unequal variance SDT mode." In S. M. Fleming & C. D. Frith (Eds.), The Cognitive Neuroscience of Metacognition (pp.25-66). Published by Springer.

If you utilize these functions, please cite the aforementioned papers and scripts upon which they are based.

The same input-format as M & L is expected. As per the original code:

% INPUTS
%
% * nR_S1, nR_S2
% these are vectors containing the total number of responses in
% each response category, conditional on presentation of S1 and S2.
%
% e.g. if nR_S1 = [100 50 20 10 5 1], then when stimulus S1 was
% presented, the subject had the following response counts:
% responded S1, rating=3 : 100 times
% responded S1, rating=2 : 50 times
% responded S1, rating=1 : 20 times
% responded S2, rating=1 : 10 times
% responded S2, rating=2 : 5 times
% responded S2, rating=3 : 1 time
%
% The ordering of response / rating counts for S2 should be the same as it
% is for S1. e.g. if nR_S2 = [3 7 8 12 27 89], then when stimulus S2 was
% presented, the subject had the following response counts:
% responded S1, rating=3 : 3 times
% responded S1, rating=2 : 7 times
% responded S1, rating=1 : 8 times
% responded S2, rating=1 : 12 times
% responded S2, rating=2 : 27 times
% responded S2, rating=3 : 89 times

The function's input should consist of counts for each of these responses separately for each stimulus type.
