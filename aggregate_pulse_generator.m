function x_hat = aggregate_pulse_generator(kj, signs, N_total)
%AGGREGATE_PULSE_GENERATOR Build the designed pulse train x_hat[k] (Eq. 9).
%
% PURPOSE
%   Assembles the sparse, alternating-sign "designed pulse train"
%   x_hat[k] from a set of pulse arrival sample indices kj and their
%   alternating signs (-1)^j, exactly as defined in Eq. (9). This is
%   the "aggregate" pulse train referred to by the modulation's name
%   (Aggregate Spread Pulse Modulation): a single sparse train
%   aggregating all the individual pulse arrivals of one transmitter's
%   payload before pulse shaping (Eq. 11-12).
%
% INPUTS
%   kj      - 1xN_sym vector of pulse arrival sample indices (0-based,
%             as produced by aspm_encoder.m / Eq. 10).
%   signs   - 1xN_sym vector of +-1 alternating signs (Eq. 9).
%   N_total - total length (in samples) of the output train x_hat.
%
% OUTPUTS
%   x_hat - 1xN_total row vector, zero everywhere except x_hat(kj+1)=signs
%           (MATLAB is 1-indexed; kj is the paper's 0-based sample index).
%
% THEORY / PAPER REFERENCE
%   Eq. (9): x_hat[k] = sum_j [[k=kj]] * (-1)^j, where [[.]] is the
%   Iverson bracket (Sec. II; D. E. Knuth, "Two notes on notation," [8]).
%   The alternating sign ensures x_hat[k] is zero-mean, eliminating a
%   DC bias in the modulating signal (Sec. II, text following Eq. 9).
%
% EQUATION REFERENCES
%   Eq. (9).

if any(kj < 0) || any(kj >= N_total)
    error('aggregate_pulse_generator:range', ...
        'All pulse indices kj must satisfy 0 <= kj < N_total.');
end

x_hat = zeros(1, N_total);
x_hat(kj + 1) = signs;  % Eq. (9), converted to MATLAB 1-based indexing

end
