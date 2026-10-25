function h = helper_functions()
%HELPER_FUNCTIONS Small shared utilities used across the M-ASPM project.
%
% PURPOSE
%   Returns a struct 'h' of function handles for small pieces of math
%   that recur throughout the project (Iverson bracket, dB conversions,
%   binomial coefficients in log-space, harmonic/arithmetic means, and
%   a generic bisection-based inverse-BER solver). Bundling them as a
%   struct of handles keeps every one-off utility out of its own file
%   while still giving each a clear, documented purpose.
%
% INPUTS
%   (none)
%
% OUTPUTS
%   h - struct of function handles, see fields below.
%
% THEORY / PAPER REFERENCE
%   h.iverson    : Iverson bracket [[.]], used in Eq. (9), (11), (13);
%                  see D. E. Knuth, "Two notes on notation," [8].
%   h.dB / h.idB : standard dB <-> linear conversions used throughout.
%   h.log_nchoosek : log of binomial coefficient (Eq. (17)), computed
%                  via gammaln to avoid overflow for large M.
%   h.harmmean_local, h.arithmean_local : harmonic/arithmetic means,
%                  used in Eqs. (43)-(52).
%
% EQUATION REFERENCES
%   Eq. (9), (11), (13), (17), (43)-(52).

h = struct();

% Iverson bracket: [[condition]] = 1 if true, 0 if false. Eqs. (9),(11),(13).
h.iverson = @(cond) double(cond);

% dB <-> linear power ratio conversions.
h.dB  = @(x) 10*log10(x);
h.idB = @(x) 10.^(x/10);

% log( C(n,k) ) via gammaln, robust for large n (used in Eq. (17)).
h.log_nchoosek = @(n,k) gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1);

% Harmonic mean of a vector (Eqs. (43),(44),(51)); avoids requiring the
% Statistics and Machine Learning Toolbox's harmmean().
h.harmmean_local = @(x) numel(x) ./ sum(1./x);

% Arithmetic mean (Eqs. (50),(52)).
h.arithmean_local = @(x) mean(x);

% Generic monotone bisection solver: find x in [lo,hi] such that
% f(x) = target, where f is assumed monotone (increasing or
% decreasing) on [lo,hi]. Used to invert the BER formula, Eq. (17)-
% (18), and to solve nonlinear partition equations, Eqs. (36),(39).
h.bisection_solve = @(f, target, lo, hi, tol) local_bisection(f, target, lo, hi, tol);

end

% ------------------------------------------------------------------
function x = local_bisection(f, target, lo, hi, tol)
% LOCAL_BISECTION Robust bisection root finder for g(x) = f(x)-target.
%   Handles both increasing and decreasing monotone f automatically.
    if nargin < 5 || isempty(tol)
        tol = 1e-10;
    end
    glo = f(lo) - target;
    ghi = f(hi) - target;
    if glo == 0, x = lo; return; end
    if ghi == 0, x = hi; return; end
    if sign(glo) == sign(ghi)
        error('helper_functions:bisection_solve', ...
            'f(lo) and f(hi) must straddle the target value.');
    end
    for iter = 1:200
        mid = 0.5*(lo+hi);
        gmid = f(mid) - target;
        if gmid == 0 || (hi-lo) < tol*max(1,abs(mid))
            x = mid;
            return;
        end
        if sign(gmid) == sign(glo)
            lo = mid; glo = gmid;
        else
            hi = mid; ghi = gmid;
        end
    end
    x = 0.5*(lo+hi);
end
