function mat = to_peak_matrix(peaks, num_cols)

validateattributes( peaks, {'double'}, {'column', 'integer'} ...
  , mfilename, 'peaks' );

mat = false( numel(peaks), num_cols );
for i = 1:numel(peaks)
  if ( peaks(i) ~= 0 )
    mat(i, peaks(i):end) = true;
  end
end

end