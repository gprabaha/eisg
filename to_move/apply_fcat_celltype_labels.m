function spike_labels = apply_fcat_celltype_labels(spike_labels, ct_labels, spike_mask)

if ( nargin < 3 )
  spike_mask = find( spike_labels, {'valid-unit', 'maybe-valid-unit'} );
end

[uuid_I, uuids] = findall( spike_labels, 'uuid', spike_mask );
match_I = bfw.find_combinations( ct_labels, uuids );

for i = 1:numel(uuid_I)
  mi = match_I{i};
  if ( ~isempty(mi) )
    assert( numel(mi) == 1, 'Expected 1 or 0 matches for unit id; got %d', numel(mi) );
    ct_label = cellstr( ct_labels, 'cell-type', mi );
    addsetcat( spike_labels, 'cell-type', ct_label, uuid_I{i} );
  end
end

end