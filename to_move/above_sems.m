function tf = above_sems(psth, psth_I, means, sigma, n)

validateattributes( means, {'double'}, {'column'}, mfilename, 'means' );
validateattributes( sigma, {'double'}, {'column'}, mfilename, 'sems' );
assert( numel(means) == numel(sigma) && numel(means) );
assert( size(psth, 1) == numel(psth_I) );

tf = false( size(psth) );

for i = 1:size(psth, 1)
  ind = psth_I{i};
  if ( ~isempty(ind) )
    assert( numel(ind) == 1, 'Expected 0 or 1 elements in index.' );
    tf(i, :) = psth(i, :) > means(ind, :) + sigma(ind, :) .* n;
  end
end

end