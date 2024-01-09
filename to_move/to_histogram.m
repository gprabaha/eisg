function ls = to_histogram(peak_mat, each_I)

if ( nargin == 1 )
  each_I = { 1:size(peak_mat, 1) };
end

ls = zeros( numel(each_I), size(peak_mat, 2) );

for i = 1:numel(each_I)
  ei = each_I{i};
  pm = double( peak_mat(ei, :) );
  sr = sum( pm, 1 );
  
  empties = all( pm == 0, 2 );
  num_non_empties = sum( ~empties );
  
  if ( num_non_empties > 0 )
    ls(i, :) = sr / num_non_empties;
  else
    ls(i, :) = mean( pm, 1 );
  end
end

end