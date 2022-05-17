%%

data_p = '/Users/prabaha/repositories/eisg/processed_data';

sorted = shared_utils.io.fload( fullfile(data_p, 'sorted_neural_data_social_gaze.mat') );
[spike_ts, spike_labels] = linearize_sorted( sorted );

unit_mask = find( spike_labels, {'valid-unit', 'maybe-valid-unit'} );

[unit_wfs, wf_labels] = extract_wfs(sorted, spike_labels, unit_mask);

[unit_p2vs, p2v_labels] = extract_p2vs(sorted, spike_labels, unit_mask);

%%

valid_unit_features = eisg.util.extract_valid_unit_features_from_sorted(sorted, params);
%% Classify celltypes based on peak-to-valley by region

regions = unique( {sorted_neural_data.region} );

n_cluster_range = [2 7];

% 'b' stands for 'broad, 'm' for 'medium', and 'n' for 'narrow' spiking
% units
[unit_celltypes, celltype_labels] = classify_regional_wfs_by_p2v

celltype_labels = eisg.celltype_class.do_xiang_multi_cluster_kmeans_by_region( valid_unit_features, regions, n_cluster_range );

label_filename = 'celltype_labels_classified_by_region.mat';
label_path = fullfile( params.processed_data_folder, label_filename );
save(label_path, 'celltype_labels');

%% Classify celltypes based on peak-to-valley across PFC and separately for BLA


%% Classify celltypes based on entire WF shape by regions


%% Classify celltypes based on entire WF across PFC and separately for BLA



%%

spike_labels = apply_cell_type_labels(...
  spike_labels, ct_labels.label_mat, ct_labels.label_mat_cols );