function [bits_hat, mj_hat] = aspm_decoder(delta_k_hat, M, n_spacing)
%ASPM_DECODER Recover bits from detected pulse-position offsets.
%
% PURPOSE
%   Inverts the encoding performed by aspm_encoder.m / pulse_position_
%   mapper.m: converts detected offsets Delta_k_hat[j] (from
%   pulse_detector.m) back into M-ary symbol indices and then into the
%   original information bits.
%
% INPUTS
%   delta_k_hat - 1xN_sym vector of estimated pulse-position offsets.
%   M           - M-ary alphabet size.
%   n_spacing   - pulse-position spacing "n" used by the encoder.
%
% OUTPUTS
%   bits_hat - 1x(N_sym*log2(M)) recovered bit sequence.
%   mj_hat   - 1xN_sym recovered symbol indices (1..M).
%
% THEORY / PAPER REFERENCE
%   Inverse of Eq. (10)/(27) via pulse_position_mapper.m's inverse
%   map, followed by converting each symbol index (1..M) back to its
%   log2(M)-bit MSB-first binary representation (inverse of the
%   encoding convention used in aspm_encoder.m).
%
% EQUATION REFERENCES
%   Inverse of Eq. (9)-(10), (27).

ppm = pulse_position_mapper();
mj_hat = ppm.inverse(delta_k_hat, n_spacing);
mj_hat = min(max(mj_hat, 1), M); % clip against detection noise pushing index out of [1,M]

bits_per_symbol = round(log2(M));
symbol_values = mj_hat - 1; % back to 0..M-1

N_sym = numel(mj_hat);
bits_matrix = zeros(bits_per_symbol, N_sym);
for j = 1:N_sym
    b = de2bi_msb(symbol_values(j), bits_per_symbol);
    bits_matrix(:, j) = b;
end

bits_hat = reshape(bits_matrix, 1, []);

end

% -------------------------------------------------------------------
function b = de2bi_msb(val, n_bits)
    % MSB-first binary representation of val using n_bits bits
    % (matches aspm_encoder.m's weights = 2.^(n_bits-1:-1:0) convention).
    b = zeros(n_bits, 1);
    for i = 1:n_bits
        shift = n_bits - i;
        b(i) = bitshift(val, -shift);
        val = val - b(i)*2^shift;
    end
end
