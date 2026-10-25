function itl = block_interleaver(nblocks, blocklen)
%BLOCK_INTERLEAVER Full-depth row/column block interleaver.
%
% NOVELTY-SUPPORT UTILITY (see README section 9).
%   An M-ary symbol error flips up to log2(M) *adjacent* coded bits at once,
%   so a short block code such as Hamming(7,4) placed directly on the
%   M-ASPM bit stream sees these errors as bursts far exceeding its t = 1
%   correction power and gains nothing. Interleaving the coded bits over
%   many codewords before transmission (and de-interleaving at the
%   receiver) spreads each burst across distinct codewords, restoring the
%   approximately-independent-bit-error condition assumed by
%   coded_ber_hamming.m and by classical block-code analysis.
%
% CONSTRUCTION
%   fwd: the length-(nblocks*blocklen) coded stream is written row-wise into
%        an nblocks x blocklen array (row = one codeword) and read out
%        column-wise, so consecutive transmitted bits belong to consecutive
%        codewords. A burst of up to nblocks consecutive channel-bit errors
%        then hits at most one bit per codeword.
%   inv: exact inverse permutation.
%
% USAGE
%   itl = block_interleaver(nblocks, hc.n);
%   tx  = itl.fwd(coded_bits);          % numel(coded_bits) == itl.len
%   coded_hat = itl.inv(rx_bits);

itl = struct();
itl.nblocks  = nblocks;
itl.blocklen = blocklen;
itl.len      = nblocks * blocklen;
itl.fwd = @(x) reshape(reshape(x(:).', blocklen, nblocks).', 1, []);
itl.inv = @(x) reshape(reshape(x(:).', nblocks, blocklen).', 1, []);

end
