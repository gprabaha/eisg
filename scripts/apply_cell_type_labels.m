function l = apply_cell_type_labels(l, ct_label_mat, ct_label_cols)

assert( size(ct_label_mat, 2) == numel(ct_label_cols) ...
  , 'Label columns do not correspond to label matrix.' );

need_cols = {'uuid', 'celltype_label'};
[tf, inds] = ismember( need_cols, ct_label_cols );
assert( all(tf), 'Missing label columns: %s', strjoin(need_cols(~tf), ' | ') );

uuid_search = cellfun( @(x) sprintf('uuid-%s', x) ...
  , cellstr(ct_label_mat(:, inds(1))), 'un', 0 );
uuid_I = cellfun( @(x) find(l, x), uuid_search, 'un', 0 );

for j = 1:numel(uuid_I)
  ct_lab = char(ct_label_mat(j, inds(2)));
  if ( ~strcmp(ct_lab, '<undefined>') )
    addsetcat( l, 'cell-type', ct_lab, uuid_I{j} );
  end
end

prune( l );

end