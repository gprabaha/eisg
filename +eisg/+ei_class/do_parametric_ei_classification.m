function do_parametric_ei_classification(sorted_neural_data, varargin)

defaults = eisg.util.make_analysis_params_struct();
defaults.feature_list = {
    'peak_to_valley', ...
    'halfwidth', ...
    'peak_trough_ratio', ...
    'recovery_slope', ...
    'repolarization_slope'
};
defaults.ei_data_filename_prefix = 'ei_classification_social_gaze';
defaults.ei_classification_method = 'kmeans';
defaults.n_ei_clusters = 2:3;
defaults.features_to_classify_by = {1, [1 2], [1 5], [1 2 3], [1 2 5], [1 2 3 4 5]};
defaults.use_calculated_features = true;
defaults.save_ei_mat = true;

params = shared_utils.general.parsestruct( defaults, varargin );

feature_list = params.feature_list;
n_ei_clusters = params.n_ei_clusters;
ei_classification_method = params.ei_classification_method;
features_to_classify_by = params.features_to_classify_by;
ei_data_filename_prefix = params.ei_data_filename_prefix;
save_ei_mat = params.save_ei_mat;

valid_unit_features = eisg.util.extract_valid_unit_features_from_sorted(sorted_neural_data, params);

region_col_in_labels = strcmp(valid_unit_features.label_categories, 'region');

ei_label_mat_cols = valid_unit_features.label_categories;
ei_label_mat_cols{end+1} = 'ei_label';
ei_label_mat_cols{end+1} = 'peak_to_valley';

valid_unit_labels = valid_unit_features.label_matrix;
regions = unique( valid_unit_labels(:, region_col_in_labels) );

progress_counter = 0;
total_steps = length(regions)*numel(features_to_classify_by)*numel(n_ei_clusters);
for region_ind = 1:length(regions)
    clc;
    disp('Parametric ei classification progress:');
    eisg.util.draw_progress_bar(progress_counter, total_steps, params.num_ticks_in_progress_bar);
    region = regions(region_ind);
    regional_unit_inds =  valid_unit_labels(:, region_col_in_labels) == region;
    regional_unit_labels = valid_unit_labels(regional_unit_inds, :);
    regional_unit_features = valid_unit_features.feature_matrix(regional_unit_inds, :);
    for feature_set_ind = 1:numel(features_to_classify_by)
        for n_clusters = n_ei_clusters
            ei_labels = struct();
            feature_cols = features_to_classify_by{feature_set_ind};
            feature_mat_for_classification = regional_unit_features(:,feature_cols);
            p2v = categorical(feature_mat_for_classification(:, 1));
            extracted_ei_labels = eisg.ei_class.generate_ei_labels_from_feature_mat(feature_mat_for_classification, n_clusters);
            ei_label_mat = [regional_unit_labels extracted_ei_labels p2v];
            ei_labels.label_mat = ei_label_mat;
            ei_labels.label_mat_cols = ei_label_mat_cols;
            ei_labels.classification_features = feature_list(feature_cols);
            ei_labels.n_clusters = n_clusters;
            if save_ei_mat
                feature_str = num2str(feature_cols);
                feature_str = strrep(feature_str, ' ', '_');
                filename = sprintf('%s-%s-%d_clusters-%s_features-%s_classification.mat', char(region), ei_data_filename_prefix, n_clusters, feature_str, ei_classification_method);
                destination_folder = fullfile(params.processed_data_folder, 'parametric_ei_labels/');
                if ~exist(destination_folder, 'dir')
                    mkdir(destination_folder);
                end
                filepath = fullfile(destination_folder, filename);
                save(filepath, 'ei_labels');
            end
            progress_counter = progress_counter + 1;
        end
    end
end
clc;
disp('Parametric ei classification progress:');
eisg.util.draw_progress_bar(progress_counter, total_steps, params.num_ticks_in_progress_bar);

end