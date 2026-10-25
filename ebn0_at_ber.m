function ebn0_dB = ebn0_at_ber(ebn0_dB_grid, ber_curve, target_ber)
%EBN0_AT_BER Eb/N0 (dB) at which a BER curve crosses a target BER, by
%   linear interpolation in (Eb/N0 dB, log10 BER). Returns NaN when the
%   target lies outside the sampled curve.
%
% Small helper for the novelty Hamming coding-gain reporting
% (plot_fig_hamming.m, main.m Part 2b). Not part of the paper.

ebn0_dB_grid = ebn0_dB_grid(:);
ber_curve    = ber_curve(:);

ok = isfinite(ber_curve) & ber_curve > 0;
ebn0_dB_grid = ebn0_dB_grid(ok);
ber_curve    = ber_curve(ok);

if isempty(ber_curve) || target_ber < min(ber_curve) || target_ber > max(ber_curve)
    ebn0_dB = NaN;
    return;
end

[ber_sorted, idx] = sort(ber_curve, 'descend');
x_sorted = ebn0_dB_grid(idx);
lb = log10(ber_sorted);

[lb_u, iu] = unique(lb, 'stable');
x_u = x_sorted(iu);

ebn0_dB = interp1(lb_u, x_u, log10(target_ber), 'linear');

end
