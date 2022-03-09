function params = make_analysis_params_struct()
 
params = struct();
 
% Common, and PSTH extraction
params.processed_data_folder = 'processed_data';
params.include_unsure_units = false;
params.sorted_data_filename = [];
params.trial_data_filename = [];
params.ei_label_filename = [];
 
% ROI-s
params.face_roi = {'face', 'eyes_nf'};
params.eye_roi = {'eyes_nf'};
params.ne_face_roi = {'face'};
params.obj_roi = {'left_nonsocial_object', 'right_nonsocial_object'};

% PSTH plotting and analysis
params.plot_folder = 'plots';
 
params.region = 'bla';
params.tasktype = 'free_viewing';
params.population_psth_method = 'mean';

% E-I classification
params.ei_labels_folder = 'parametric_ei_labels';
params.ei_data_filename_prefix = [];
params.ei_classification_method = 'kmeans';
params.feature_list = {
    'peak_to_valley', ...
    'halfwidth', ...
    'peak_trough_ratio', ...
    'recovery_slope', ...
    'repolarization_slope'
};
params.n_ei_clusters = 2;
params.use_calculated_features = true;
params.features_to_classify_by = {[1, 2]};
params.save_ei_mat = true;
params.default_ei_features = [1, 2];
params.default_ei_clusters = 2;
 
% Progress bar
params.num_ticks_in_progress_bar = 75;
 
end
