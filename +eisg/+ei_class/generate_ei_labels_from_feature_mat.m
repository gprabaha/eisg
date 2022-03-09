function ei_labels = generate_ei_labels_from_feature_mat(feature_mat, n_clusters)

[idx, ~] = kmeans(feature_mat, n_clusters);
max_val = -inf;
min_val = inf;
max_cluster = nan;
min_cluster = nan;
p2v_index = 1;
ei_labels = repmat({'o'}, [length(idx), 1]);
for cluster_index = 1:n_clusters
   mean_p2v = nanmean( feature_mat(idx==cluster_index, p2v_index) );
   if mean_p2v > max_val
       max_val = mean_p2v;
       max_cluster = cluster_index;
   end
   if mean_p2v < min_val
       min_val = mean_p2v;
       min_cluster = cluster_index;
   end
end
ei_labels(idx==max_cluster) = {'e'};
ei_labels(idx==min_cluster) = {'i'};

end