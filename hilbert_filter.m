function hf = hilbert_filter()
%HILBERT_FILTER Discrete Hilbert transform used to build the PSF's
%   imaginary (quadrature) component from its real component.
%
% PURPOSE
%   Returns a struct 'hf' of function handles implementing the
%   relation h[k] = H{g[k]} used in Eq. (11): the imaginary part of a
%   pulse-shaping filter (PSF) impulse response is the discrete
%   Hilbert transform of its real part, so that the resulting complex
%   pulse occupies a single (analytic-signal) sideband once mixed to a
%   carrier (Sec. II, text following Eq. (11); Sec. II, "Since h_i[k]
%   is the Hilbert transform of g_i[k], this waveform will occupy only
%   a single sideband...[9]").
%
% INPUTS
%   (none)
%
% OUTPUTS
%   hf - struct of function handles, described below.
%
% THEORY / PAPER REFERENCE
%   Eq. (11): h[k] = H{g[k]}. The paper cites R. N. Bracewell, "The
%   Fourier Transform and Its Applications," and Todoran et al.,
%   "Discrete Hilbert transform, numeric algorithms," as [9],[10] for
%   the discrete Hilbert transform itself. We implement it using
%   MATLAB's built-in hilbert() (Signal Processing Toolbox), which
%   returns the analytic signal g[k] + i*H{g[k]} via the standard
%   FFT-based method described in exactly those references; taking the
%   imaginary part gives H{g[k]} directly.
%
% EQUATION REFERENCES
%   Eq. (11).

hf = struct();
hf.transform = @(g) imag(hilbert(g(:).'));

end
