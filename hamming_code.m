function hc = hamming_code(m)
%HAMMING_CODE Systematic single-error-correcting Hamming(2^m-1, 2^m-1-m) codec.
%
% NOVELTY ADDITION (not part of Nikitin & Davidchack, IEEE Access 2023).
%   The paper analyses UNCODED noncoherent M-ASPM (Eq. 17, Figs. 3-4). This
%   module adds a classical outer forward-error-correction (FEC) layer -- a
%   binary Hamming block code -- around the paper's pulse-position bit stream,
%   so the effect of a low-complexity FEC on the M-ASPM BER / range / SIR
%   trade-off can be quantified. Hamming was chosen because it is the
%   lowest-complexity perfect code (hard-decision syndrome decoding, no
%   toolbox, negligible receiver cost) -- a natural fit for the paper's
%   "highly scalable, low-power" LPWAN target.
%
% CODE
%   n = 2^m - 1 coded bits, k = n - m information bits, minimum distance 3,
%   corrects any single bit error per block (t = 1). Standard choices:
%     m = 3  -> Hamming(7,4),   rate 0.571
%     m = 4  -> Hamming(15,11), rate 0.733
%     m = 5  -> Hamming(31,26), rate 0.839
%   Layout is the classical Hamming arrangement: parity bits occupy the
%   power-of-two positions 1,2,4,...,2^(m-1); the parity-check matrix column
%   for position j is the m-bit binary representation of j, so a nonzero
%   syndrome read as an integer IS the 1-indexed position of the single error.
%
% USAGE
%   hc = hamming_code(3);
%   coded = hc.encode(info_bits);     % numel(info_bits) must be a multiple of hc.k
%   info_hat = hc.decode(rx_bits);    % numel(rx_bits)   must be a multiple of hc.n
%
% OUTPUT STRUCT FIELDS
%   .m .n .k .rate            code parameters
%   .parity_pos .data_pos     1-indexed bit positions within a block
%   .H                        m x n parity-check matrix (GF(2))
%   .encode(bits)             info bits -> systematic codeword bits
%   .decode(bits)             received bits -> corrected info bits
%
% Requires no MATLAB toolboxes (pure GF(2) linear algebra with mod-2).

if nargin < 1 || isempty(m), m = 3; end
if m < 2 || m ~= round(m)
    error('hamming_code:m', 'm must be an integer >= 2.');
end

n = 2^m - 1;
k = n - m;
parity_pos = 2.^(0:m-1);                 % 1-indexed positions holding parity bits
data_pos   = setdiff(1:n, parity_pos);   % 1-indexed positions holding data bits (ascending)

% H(i,j) = bit (i-1) of the integer j  ->  column j is the binary code of j.
cols = 1:n;
H = rem(floor(cols ./ (2.^((0:m-1)'))), 2);   % m x n, entries in {0,1}

hc = struct();
hc.m = m; hc.n = n; hc.k = k; hc.rate = k / n;
hc.parity_pos = parity_pos;
hc.data_pos   = data_pos;
hc.H = H;
hc.encode = @(bits) local_encode(bits, n, k, m, parity_pos, data_pos, H);
hc.decode = @(bits) local_decode(bits, n, m, data_pos, H);

end

% -------------------------------------------------------------------
function coded = local_encode(bits, n, k, m, parity_pos, data_pos, H)
    bits = bits(:).';
    if mod(numel(bits), k) ~= 0
        error('hamming_code:encodeLen', ...
            'Number of information bits (%d) must be a multiple of k = %d.', numel(bits), k);
    end
    nblk = numel(bits) / k;
    D = reshape(bits, k, nblk);          % k x nblk, one info block per column

    C = zeros(n, nblk);
    C(data_pos, :) = D;
    % Choose each parity bit so its parity-check row sums to zero (mod 2).
    for i = 1:m
        C(parity_pos(i), :) = mod(H(i, data_pos) * D, 2);
    end
    coded = reshape(C, 1, []);
end

% -------------------------------------------------------------------
function info = local_decode(bits, n, m, data_pos, H)
    bits = bits(:).';
    if mod(numel(bits), n) ~= 0
        error('hamming_code:decodeLen', ...
            'Number of received bits (%d) must be a multiple of n = %d.', numel(bits), n);
    end
    nblk = numel(bits) / n;
    R = reshape(bits, n, nblk);          % n x nblk, one received block per column

    S = mod(H * R, 2);                    % m x nblk syndrome
    synd_int = (2.^(0:m-1)) * S;          % 1 x nblk, 0 (no error) .. n (error position)

    for b = 1:nblk
        e = synd_int(b);
        if e >= 1 && e <= n
            R(e, b) = 1 - R(e, b);       % correct the single detected error
        end
    end

    info = reshape(R(data_pos, :), 1, []);
end
