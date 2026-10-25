function mc = multi_user_channel()
%MULTI_USER_CHANNEL Combine multiple M-ASPM transmitters at a receiver.
%
% PURPOSE
%   Returns a struct 'mc' with a function that sums several
%   transmitters' baseband waveforms (each with its own PSF, path-loss
%   attenuation and time offset) into one composite received signal,
%   reproducing the inter-PSF mutual-interference scenario illustrated
%   in Fig. 9 (Sec. III-B).
%
% INPUTS
%   (none) -- see mc.combine(...) below.
%
% OUTPUTS
%   mc - struct with field 'combine', a function handle:
%        r = mc.combine(signals, amplitudes, delays)
%
% THEORY / PAPER REFERENCE
%   Section III-B, Fig. 9: multiple equal-IpI M-ASPM transmitters at
%   different ranges d_i, with pulse duty cycles D_i chosen so that
%   D_i/D_j = (d_i/d_j)^gamma (Eq. 25 with N_i=N_j), producing equal
%   AVERAGE RECEIVED POWER for every node (so that any residual SIR
%   difference among nodes comes only from interference statistics,
%   not from unequal desired-signal strength -- Sec. III-B text).
%   The received composite signal is simply the sum of the individual
%   (already amplitude-scaled) transmit waveforms, each convolved with
%   its own channel delay -- standard superposition for a shared,
%   linear multiple-access AWGN channel.
%
% EQUATION REFERENCES
%   Eq. (23)-(25) (relative attenuation via duty cycle / path loss).

mc = struct();
mc.combine = @(signals, amplitudes, delays) local_combine(signals, amplitudes, delays);

end

% -------------------------------------------------------------------
function r = local_combine(signals, amplitudes, delays)
    % signals    - 1xN_tx cell array of complex row-vector waveforms
    % amplitudes - 1xN_tx vector of linear amplitude scale factors
    %              (e.g. sqrt(received power) relative to a reference)
    % delays     - 1xN_tx vector of integer sample delays (>=0)

    N_tx = numel(signals);
    max_len = 0;
    for i = 1:N_tx
        max_len = max(max_len, delays(i) + numel(signals{i}));
    end

    r = zeros(1, max_len);
    for i = 1:N_tx
        idx = delays(i) + (1:numel(signals{i}));
        r(idx) = r(idx) + amplitudes(i) .* signals{i};
    end
end
