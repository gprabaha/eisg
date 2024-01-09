function [rs_tables, rs_labels] = ranksums(data, labels, each, as, bs, label_category)

if ( nargin < 6 )
  label_category = [];
end

assert( numel(as) == numel(bs) );

rs_tables = {};
rs_labels = fcat();

for i = 1:numel(as)
  a = as{i};
  b = bs{i};
  
  rs_outs = dsp3.ranksum( data, labels, each, a, b );
  if ( ~isempty(label_category) )
    al = strjoin( cellstr(a), ' | ' );
    bl = strjoin( cellstr(b), ' | ' );
    setcat( rs_outs.rs_labels, label_category, sprintf('%s_%s', al, bl) );
  end
  
  rs_tables = [ rs_tables; rs_outs.rs_tables ];
  append( rs_labels, rs_outs.rs_labels );
end

end