function l = apply_face_obj_sig_labels(l, sig_cells)

region_col = 1;
uuid_col = 2;

uuid_search = sig_cells(:, uuid_col);
uuid_I = cellfun( @(x) find(l, x), uuid_search, 'un', 0 );

for j = 1:numel(uuid_I)
  ct_lab = 'sig-f-o';
  addsetcat( l, 'face-vs-obj', ct_lab, uuid_I{j} );
end

prune( l );

end