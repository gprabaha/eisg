function celltype_labels = classify_regional_wfs_by_p2v( valid_unit_features, regions, n_cluster_range )

if nargin < 3
    minNumOfCluster = 2;
    maxNumOfCluster = 9;
else
    minNumOfCluster = n_cluster_range(1);
    maxNumOfCluster = n_cluster_range(2);
end
offset = minNumOfCluster - 1;

feature_mat = valid_unit_features.feature_matrix;
feature_list = valid_unit_features.feature_list;
label_matrix = valid_unit_features.label_matrix;
label_categories = valid_unit_features.label_categories;

region_col_in_label_mat = strcmp( label_categories, 'region' );
all_regions = unique( label_matrix(:,region_col_in_label_mat) );

regions_to_eliminate = setdiff( cellstr( all_regions' ), regions );
for region = regions_to_eliminate
    region_inds = label_matrix(:,region_col_in_label_mat) == region;
    feature_mat(region_inds,:) = [];
    label_matrix(region_inds, :) = [];
end
celltype_label_mat = label_matrix;

celltype_labels = struct();
celltype_label_mat_cols = valid_unit_features.label_categories;
celltype_label_mat_cols{end+1} = 'celltype_label';
celltype_label_mat_cols{end+1} = 'peak_to_valley_mu_s';

n_clusters = [];
for region = regions
    disp(['running for region:' char(region)]);
    p2v_col_in_feature_mat = strcmp( feature_list, 'peak_to_valley' );
    region_col_in_label_mat = strcmp( label_categories, 'region' );
    region_inds = label_matrix(:, region_col_in_label_mat) == region;
    regional_p2vs = feature_mat(region_inds, p2v_col_in_feature_mat);
    p2v_mu_s = categorical( round ( regional_p2vs * 1e6 ) );
    celltype_label_mat(region_inds, 6) = p2v_mu_s;
    region_inds = find( region_inds );
    % Xiang's part
    datanonnorm4Cluster = standardization(regional_p2vs);
    [clustFilt,labels,percElements,centroids] = kMeansClustering(datanonnorm4Cluster,minNumOfCluster,maxNumOfCluster);
    [aic,bic,fom_aic,fom_bic,nclus_aic,nclus_bic] = AICBICAnalysis(datanonnorm4Cluster, labels, clustFilt, minNumOfCluster, maxNumOfCluster, offset);
    rel_nclus_ind = nclus_aic - offset;
    optimal_clustFilt = clustFilt(rel_nclus_ind, :);
    optimal_centroids = centroids(rel_nclus_ind, :);
    n_units_in_each_cluster = cellfun(@(x) numel(x), optimal_clustFilt);
    
    if nclus_aic > 2
        n_top_clusters = 3;
        [~,top_cluster_inds] = maxk(n_units_in_each_cluster, n_top_clusters);
        top_centroids = cell2mat( optimal_centroids(top_cluster_inds) );

        [~, narrow_ind] = min(top_centroids);
        narrow_cluster = top_cluster_inds(narrow_ind);
        rel_narrow_unit_inds = optimal_clustFilt{narrow_cluster};
        narrow_unit_inds = region_inds(rel_narrow_unit_inds);
        celltype_label_mat(narrow_unit_inds, 5) = 'n';

        [~, broad_ind] = max(top_centroids);
        broad_cluster = top_cluster_inds(broad_ind);
        rel_broad_unit_inds = optimal_clustFilt{broad_cluster};
        broad_unit_inds = region_inds(rel_broad_unit_inds);
        celltype_label_mat(broad_unit_inds, 5) = 'b';

        medium_ind = setdiff(1:3, [narrow_ind, broad_ind]);
        medium_cluster = top_cluster_inds(medium_ind);
        rel_medium_unit_inds = optimal_clustFilt{medium_cluster};
        medium_unit_inds = region_inds(rel_medium_unit_inds);
        celltype_label_mat(medium_unit_inds, 5) = 'm';

    else
        n_top_clusters = 2;
        [~,top_cluster_inds] = maxk(n_units_in_each_cluster, n_top_clusters);
        top_centroids = cell2mat( optimal_centroids(top_cluster_inds) );

        [~, narrow_ind] = min(top_centroids);
        narrow_cluster = top_cluster_inds(narrow_ind);
        rel_narrow_unit_inds = optimal_clustFilt{narrow_cluster};
        narrow_unit_inds = region_inds(rel_narrow_unit_inds);
        celltype_label_mat(narrow_unit_inds, 5) = 'n';

        [~, broad_ind] = max(top_centroids);
        broad_cluster = top_cluster_inds(broad_ind);
        rel_broad_unit_inds = optimal_clustFilt{broad_cluster};
        broad_unit_inds = region_inds(rel_broad_unit_inds);
        celltype_label_mat(broad_unit_inds, 5) = 'b';

    end
    n_clusters = [n_clusters nclus_aic];
end

celltype_labels.label_mat = celltype_label_mat;
celltype_labels.label_mat_cols = celltype_label_mat_cols;
celltype_labels.regions = regions;
celltype_labels.n_clusters_per_region = n_clusters;


end