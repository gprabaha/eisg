function [coh, coh_labels, f, t] = load_tf_measure(ms, mean_each)

coh_labels = fcat();
coh = cell( numel(ms), 1 );
f = cell( size(coh) );
t = cell( size(coh) );

for i = 1:numel(ms)  
  fprintf( '%d of %d\n', i, numel(ms) );
  coh_file = load( ms{i} );
  tmp_coh_labels = fcat.from( coh_file.var );
  tmp_coh = coh_file.var.coh;
  
  if ( isempty(mean_each) )
    mean_coh = tmp_coh;
    mean_labs = tmp_coh_labels;
  else
    [mean_labs, mean_I] = keepeach( tmp_coh_labels', mean_each );
    mean_coh = bfw.row_nanmean( double(tmp_coh), mean_I );
  end
  
  coh{i} = mean_coh;
  f{i} = coh_file.var.f;
  t{i} = coh_file.var.t;
  append( coh_labels, mean_labs );
end

coh = vertcat( coh{:} );
assert_ispair( coh, coh_labels );

if ( isempty(f) )
  f = [];
else
  f = f{1};
end
if ( isempty(t) )
  t = [];
else
  t = t{1};
end

end