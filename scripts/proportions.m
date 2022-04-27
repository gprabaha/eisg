function prop = proportions(tf, I)

prop = zeros( size(I) );
for i = 1:numel(I)
  prop(i) = pnz( tf(I{i}) );
end

end