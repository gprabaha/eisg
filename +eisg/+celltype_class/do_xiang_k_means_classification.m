function [unit_celltypes, celltype_labels] = do_xiang_k_means_classification(class_feature_matrix, class_feature_labels, p2v_vec, n_cluster_range, use_pca)

if nargin < 5
  use_pca = false;
end
if nargin < 4
    minNumOfCluster = 2;
    maxNumOfCluster = 9;
else
    minNumOfCluster = n_cluster_range(1);
    maxNumOfCluster = n_cluster_range(2);
end
offset = minNumOfCluster - 1;

if use_pca
  % Get the first PC
  if size(class_feature_matrix, 2) > 1
    mat_for_kmeans_classification = pca( class_feature_matrix' );
    mat_for_kmeans_classification = mat_for_kmeans_classification(:, 1:min( size( class_feature_matrix, 2 ) - 1, 3 ));
  else
    mat_for_kmeans_classification = class_feature_matrix;
  end
else
  mat_for_kmeans_classification = class_feature_matrix;
end
  
unit_celltypes = repmat({'<cell-type>'}, size(mat_for_kmeans_classification, 1), 1);
celltype_labels = class_feature_labels;

datanonnorm4Cluster = standardization(mat_for_kmeans_classification);
[clustFilt,labels,percElements,centroids] = kMeansClustering(datanonnorm4Cluster,minNumOfCluster,maxNumOfCluster);
[aic,bic,fom_aic,fom_bic,nclus_aic,nclus_bic] = AICBICAnalysis(datanonnorm4Cluster, labels, clustFilt, minNumOfCluster, maxNumOfCluster, offset);
rel_nclus_ind = nclus_aic - offset;
optimal_clustFilt = clustFilt(rel_nclus_ind, :);
for i = 1:numel(optimal_clustFilt)
  cluster_p2vs(i) = mean( p2v_vec(optimal_clustFilt{i}) );
end
optimal_centroids = centroids(rel_nclus_ind, :);
n_units_in_each_cluster = cellfun(@(x) numel(x), optimal_clustFilt);

if nclus_aic > 2
  n_top_clusters = 3;
  [~,top_cluster_inds] = maxk(n_units_in_each_cluster, n_top_clusters);
  optimal_centroids = optimal_centroids(top_cluster_inds);
  top_clusters_mean_p2vs = cluster_p2vs(top_cluster_inds);

  [~, narrow_ind] = min(top_clusters_mean_p2vs);
  narrow_cluster = top_cluster_inds(narrow_ind);
  narrow_unit_inds = optimal_clustFilt{narrow_cluster};
  unit_celltypes(narrow_unit_inds) = {'n'};

  [~, broad_ind] = max(top_clusters_mean_p2vs);
  broad_cluster = top_cluster_inds(broad_ind);
  broad_unit_inds = optimal_clustFilt{broad_cluster};
  unit_celltypes(broad_unit_inds) = {'b'};

  medium_ind = setdiff(1:3, [narrow_ind, broad_ind]);
  medium_cluster = top_cluster_inds(medium_ind);
  medium_unit_inds = optimal_clustFilt{medium_cluster};
  unit_celltypes(medium_unit_inds) = {'m'};
else
  n_top_clusters = 2;
  [~,top_cluster_inds] = maxk(n_units_in_each_cluster, n_top_clusters);
  optimal_centroids = optimal_centroids(top_cluster_inds);
  top_clusters_mean_p2vs = cluster_p2vs(top_cluster_inds);

  [~, narrow_ind] = min(top_clusters_mean_p2vs);
  narrow_cluster = top_cluster_inds(narrow_ind);
  narrow_unit_inds = optimal_clustFilt{narrow_cluster};
  unit_celltypes(narrow_unit_inds) = {'n'};

  [~, broad_ind] = max(top_clusters_mean_p2vs);
  broad_cluster = top_cluster_inds(broad_ind);
  broad_unit_inds = optimal_clustFilt{broad_cluster};
  unit_celltypes(broad_unit_inds) = {'b'};
end

end