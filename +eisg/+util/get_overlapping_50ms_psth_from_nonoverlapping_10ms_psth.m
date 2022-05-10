function overlapping_psth = get_overlapping_50ms_psth_from_nonoverlapping_10ms_psth(nonoverlapping_psth, method)

if nargin < 2
  method = 'square_wave';
end
switch method
  case 'square_wave'
    kernel = ones(1, 5);
  case 'gaussian'
    % accept points from 2 neighboring bins, and the bins in the boundary
    % are considered to be half standard deviation away
    x = linspace( -1, 1, 5 );
    sigma = 1;
    mu = 0;
    kernel = gaussmf(x, [sigma, mu]);
end
kernel = kernel/sum(kernel);
overlapping_psth = conv2(nonoverlapping_psth, kernel, 'same');

end