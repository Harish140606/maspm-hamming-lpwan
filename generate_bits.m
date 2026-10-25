function bits = generate_bits(n_bits, seed)
%GENERATE_BITS Random source bit generator for M-ASPM payloads.
%
% PURPOSE
%   Produces the random information bit sequence that feeds the
%   M-ASPM encoder (aspm_encoder.m), representing the uplink node's
%   raw payload before Eq. (9)-(10) pulse-position encoding.
%
% INPUTS
%   n_bits - number of bits to generate.
%   seed   - (optional) RNG seed for reproducibility. If omitted, the
%            current global RNG state is used.
%
% OUTPUTS
%   bits - 1xn_bits row vector of 0/1 values (double).
%
% THEORY / PAPER REFERENCE
%   The paper does not model the information source explicitly; it
%   assumes an arbitrary bit sequence is available for pulse-position
%   encoding (Sec. II, discussion preceding Eq. 9: "we can encode
%   information in the arrival times of the pulses"). A uniform i.i.d.
%   random bit source is the standard, uncontroversial choice for
%   simulating an uncoded digital communication payload.
%
% EQUATION REFERENCES
%   Feeds Eq. (9)-(10) (via aspm_encoder.m).

if nargin >= 2 && ~isempty(seed)
    rng(seed);
end

bits = randi([0 1], 1, n_bits);

end
