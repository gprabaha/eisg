function celltype_labels = do_shape_based_wf_class_by_subset(...
  unit_wfs, spike_labels, subset_inds, unit_wf_features, p2v_inds, n_cluster_range, use_pca)

if nargin < 7
  use_pca = false;
end

celltype_labels = fcat();

for i = 1:numel(subset_inds)
  [uuid_I, ~] = findall( spike_labels(subset_inds{i}), 'uuid' );
  uuid_I = cellfun( @(x) x(1), uuid_I );
  subset_wfs = unit_wfs( subset_inds{i} );
  classification_labels = spike_labels(uuid_I);
  
  p2v_vec = unit_wf_features(p2v_inds{i});
  
  [regional_unit_celltypes, regional_celltype_labels] = eisg.celltype_class.do_xiang_k_means_classification(subset_wfs, classification_labels, p2v_vec, n_cluster_range, use_pca);
  addsetcat( regional_celltype_labels, 'cell-type', cellstr(regional_unit_celltypes) );
  append( celltype_labels, regional_celltype_labels );
  
end

end