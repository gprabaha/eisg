function [aligned_lfp, ok_t] = align_lfp_window(channel, t0, bin_width, fs)

ib_ind = @(i) i >= 0 & i < numel(channel);

lfp_bin_size = floor( bin_width * fs );
trunc_t0 = floor( t0 ./ (1 / fs) ) * (1 / fs);
lfp_i0 = floor( trunc_t0 * fs );
lfp_i1 = lfp_i0 + lfp_bin_size;

aligned_lfp = nan( numel(t0), lfp_bin_size );

ok_t = find( ~isnan(t0) & ib_ind(lfp_i0) & ib_ind(lfp_i1) );
for i = 1:numel(ok_t)
  i0 = lfp_i0(ok_t(i)) + 1;
  i1 = lfp_i1(ok_t(i)) + 1;
  aligned_lfp(ok_t(i), :) = channel(i0:i1-1);
end

end