function [celltype_labels, full_feature_mat] = do_feature_based_wf_class_for_each_subset(...
  unit_wf_features, feature_labels, subset_feature_inds, num_features, p2v_inds, n_cluster_range, use_pca)

if nargin < 7
  use_pca = false;
end

celltype_labels = fcat();
full_feature_mat = [];

for i = 1:numel(subset_feature_inds)
  regional_labels = fcat.from( feature_labels( subset_feature_inds{i}, : ), getcats(feature_labels) );
  rmcat(regional_labels, 'wf_feature');
  [uuid_I, uuids] = findall( regional_labels, 'uuid' );
  num_units = numel( uuids );
  uuid_I = cellfun( @(x) x(1), uuid_I );
  regional_features = unit_wf_features( subset_feature_inds{i} );
  
  feature_mat = reshape( regional_features, [num_units, num_features] );
  full_feature_mat = [full_feature_mat; feature_mat];
  classification_labels = regional_labels(uuid_I);
  p2v_vec = unit_wf_features(p2v_inds{i});
  
  [regional_unit_celltypes, regional_celltype_labels] = eisg.celltype_class.do_xiang_k_means_classification(feature_mat, classification_labels, p2v_vec, n_cluster_range, use_pca);
  addsetcat( regional_celltype_labels, 'cell-type', cellstr(regional_unit_celltypes) );
  append( celltype_labels, regional_celltype_labels );
  
end

end