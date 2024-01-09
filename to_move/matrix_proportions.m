function prop = matrix_proportions(tf, I)

prop = zeros( numel(I), size(tf, 2) );
for i = 1:numel(I)
  prop(i, :) = sum( tf(I{i}, :) ) ./ numel( I{i} );
end

end