function inds = find_first(m)

inds = zeros( size(m, 1), 1 );
for i = 1:size(m, 1)
  fi = find( m(i, :), 1 );
  if ( ~isempty(fi) )
    inds(i) = fi;
  end
end

end