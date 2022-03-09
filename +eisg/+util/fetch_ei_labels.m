function ei_labels = fetch_ei_labels(region, varargin)

defaults = eisg.util.make_analysis_params_struct();
defaults.ei_data_filename_prefix = 'ei_classification_social_gaze';
defaults.save_ei_mat = true;

params = shared_utils.general.parsestruct( defaults, varargin );

ei_classification_method = params.ei_classification_method;
default_ei_features = params.default_ei_features;
ei_data_filename_prefix = params.ei_data_filename_prefix;
n_clusters = params.default_ei_clusters;

feature_str = num2str(default_ei_features);
feature_str = strrep(feature_str, ' ', '_');
filename = sprintf('%s-%s-%d_clusters-%s_features-%s_classification.mat', char(region), ei_data_filename_prefix, n_clusters, feature_str, ei_classification_method);
destination_folder = fullfile(params.processed_data_folder, params.ei_labels_folder);
filepath = fullfile(destination_folder, filename);
ei_labels = load( filepath );
ei_labels = ei_labels.ei_labels;

end