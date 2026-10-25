function pd = pulse_detector()
%PULSE_DETECTOR Peak (arrival-time) detector for the M-ASPM matched-filter output.
%
% PURPOSE
%   Returns a struct 'pd' with a function that scans the noncoherent
%   detection statistic y2[k] (from matched_filter_receiver.m) inside
%   each expected Np-sample window and returns the sample offset of
%   the largest peak, i.e. the estimated Delta_k[mj] for that pulse.
%
% INPUTS
%   (none) -- see pd.find_peaks(...) below.
%
% OUTPUTS
%   pd - struct with field 'find_peaks', a function handle:
%        delta_k_hat = pd.find_peaks(y2, N_sym, Np, deltaN, filter_delay)
%
% THEORY / PAPER REFERENCE
%   Section II, text following Eq. (12): "after synchronization we
%   would need to obtain only M=8 samples per pulse...Out of each 8
%   samples of y_nc[k], the position of the sample with the largest
%   magnitude will correspond to the position of the respective pulse
%   in the designed train." This function implements exactly that
%   maximum-magnitude search, generalized to any M via the pulse-
%   position window of width Np (since 0 <= Delta_k[m] < Np, Eq. 10).
%
%   filter_delay accounts for the group delay introduced by 'full'
%   convolution in matched_filter_receiver.m (the matched filter has
%   length L, so a pulse transmitted at 0-based sample kj produces its
%   matched-filter peak at 0-based output sample kj + (L-1)); this
%   project verifies the exact offset by a noiseless round-trip test
%   (see README.md, "Numerical Validation").
%
% EQUATION REFERENCES
%   Section II text following Eq. (12); used to invert Eq. (10).

pd = struct();
pd.find_peaks = @(y2, N_sym, Np, deltaN, filter_delay) ...
    local_find_peaks(y2, N_sym, Np, deltaN, filter_delay);

end

% -------------------------------------------------------------------
function delta_k_hat = local_find_peaks(y2, N_sym, Np, deltaN, filter_delay)
    delta_k_hat = zeros(1, N_sym);
    for j = 0:(N_sym-1)
        win_start0 = j*Np + deltaN + filter_delay;      % 0-based
        idx0 = win_start0 : (win_start0 + Np - 1);       % 0-based candidate samples
        idx1 = idx0 + 1;                                  % MATLAB 1-based indices into y2
        idx1 = idx1(idx1 >= 1 & idx1 <= numel(y2));
        if isempty(idx1)
            delta_k_hat(j+1) = NaN;
            continue;
        end
        [~, rel] = max(y2(idx1));
        peak_pos0 = idx0(1) + (rel-1);                    % 0-based absolute sample
        delta_k_hat(j+1) = peak_pos0 - win_start0;         % offset within window (0..Np-1)
    end
end
