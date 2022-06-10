function [unit_celltypes, celltype_labels] = classify_regional_wfs_by_features(class_feature_matrix, class_feature_labels, n_cluster_range)

if nargin < 3
    minNumOfCluster = 2;
    maxNumOfCluster = 9;
else
    minNumOfCluster = n_cluster_range(1);
    maxNumOfCluster = n_cluster_range(2);
end
offset = minNumOfCluster - 1;

unit_celltypes = repmat(' ', size(class_feature_matrix, 1), 1);
celltype_labels = class_feature_labels;

% Find all regions from class_feature_labels

for region = regions'
  region_inds = find( class_feature_labels, region );
  features_to_classify_by = class_feature_matrix(region_inds, :);
  % Xiang's part
  datanonnorm4Cluster = standardization(features_to_classify_by);
  [clustFilt,labels,percElements,centroids] = kMeansClustering(datanonnorm4Cluster,minNumOfCluster,maxNumOfCluster);
  [aic,bic,fom_aic,fom_bic,nclus_aic,nclus_bic] = AICBICAnalysis(datanonnorm4Cluster, labels, clustFilt, minNumOfCluster, maxNumOfCluster, offset);
  rel_nclus_ind = nclus_aic - offset;
  optimal_clustFilt = clustFilt(rel_nclus_ind, :);
  optimal_centroids = centroids(rel_nclus_ind, :);
  n_units_in_each_cluster = cellfun(@(x) numel(x), optimal_clustFilt);
  region_inds_in_full_feature_mat = find(region_inds);
  
  if nclus_aic > 2
    n_top_clusters = 3;
    [~,top_cluster_inds] = maxk(n_units_in_each_cluster, n_top_clusters);
    top_centroids = cell2mat( optimal_centroids(top_cluster_inds) );

    [~, narrow_ind] = min(top_centroids);
    narrow_cluster = top_cluster_inds(narrow_ind);
    rel_narrow_unit_inds = optimal_clustFilt{narrow_cluster};
    narrow_unit_inds = region_inds_in_full_feature_mat(rel_narrow_unit_inds);
    unit_celltypes(narrow_unit_inds) = 'n';

    [~, broad_ind] = max(top_centroids);
    broad_cluster = top_cluster_inds(broad_ind);
    rel_broad_unit_inds = optimal_clustFilt{broad_cluster};
    broad_unit_inds = region_inds_in_full_feature_mat(rel_broad_unit_inds);
    unit_celltypes(broad_unit_inds) = 'b';

    medium_ind = setdiff(1:3, [narrow_ind, broad_ind]);
    medium_cluster = top_cluster_inds(medium_ind);
    rel_medium_unit_inds = optimal_clustFilt{medium_cluster};
    medium_unit_inds = region_inds_in_full_feature_mat(rel_medium_unit_inds);
    unit_celltypes(medium_unit_inds) = 'm';
  else
    n_top_clusters = 2;
    [~,top_cluster_inds] = maxk(n_units_in_each_cluster, n_top_clusters);
    top_centroids = cell2mat( optimal_centroids(top_cluster_inds) );

    [~, narrow_ind] = min(top_centroids);
    narrow_cluster = top_cluster_inds(narrow_ind);
    rel_narrow_unit_inds = optimal_clustFilt{narrow_cluster};
    narrow_unit_inds = region_inds_in_full_feature_mat(rel_narrow_unit_inds);
    unit_celltypes(narrow_unit_inds) = 'n';

    [~, broad_ind] = max(top_centroids);
    broad_cluster = top_cluster_inds(broad_ind);
    rel_broad_unit_inds = optimal_clustFilt{broad_cluster};
    broad_unit_inds = region_inds_in_full_feature_mat(rel_broad_unit_inds);
    unit_celltypes(broad_unit_inds) = 'b';
  end
end

end