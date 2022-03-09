

%% Create the sorted neural data structure from the sorter output

params = eisg.util.make_preprocess_params_struct();

params.processed_data_folder = 'processed_data';
params.include_unsure_units = true;
params.path_to_sorter_output = '/Volumes/ExtSSD/sorted_neural_data/social_gaze/';
params.sorted_neural_data_filename = 'sorted_neural_data_social_gaze.mat';
params.path_to_validity = fullfile(params.processed_data_folder, 'unit_validity_social_gaze.mat');
params.integrate_validity = true;
params.calculate_features = true;
params.feature_list = {
    'peak_to_valley', ...
    'halfwidth', ...
    'peak_trough_ratio', ...
    'recovery_slope', ...
    'repolarization_slope'
};
params.n_extra_pts_per_gap_for_interpolation = 2;
params.num_ticks_in_progress_bar = 75;

sorted_neural_data = eisg.preprocess.create_sorted_neural_data(params);
folder_path_sorted = fullfile(params.processed_data_folder);
if ~exist(folder_path_sorted, 'dir')
    mkdir(folder_path_sorted);
end
filepath_sorted = fullfile(folder_path_sorted,...
    params.sorted_neural_data_filename);
disp(['Saving sorted neural data to: ' filepath_sorted]);
save(filepath_sorted, 'sorted_neural_data', '-v7.3');
disp('Saved');

%% Parameter initiation

params = eisg.util.make_analysis_params_struct();
params.sorted_data_filename = 'sorted_neural_data_social_gaze.mat';
params.trial_data_filename = 'trial_events_social_gaze.mat';
params.ei_label_filename = 'bla-ei_classification_social_gaze-2_clusters-1__2__5_features-kmeans_classification.mat';
params.region = 'bla';
% EI Classification
params.feature_list = {
    'peak_to_valley', ...
    'halfwidth', ...
    'peak_trough_ratio', ...
    'recovery_slope', ...
    'repolarization_slope'
};
params.ei_data_filename_prefix = 'ei_classification_social_gaze';
params.ei_classification_method = 'kmeans';
params.n_ei_clusters = 2:3;
params.features_to_classify_by = {1, [1 2], [1 5], [1 2 3], [1 2 5], [1 2 3 4 5]};
params.save_ei_mat = true;

params.include_unsure_units = true;


%% Loading sorted neural data

path_to_sorted = fullfile( params.processed_data_folder, params.sorted_data_filename );
disp('Loading sorted neural data');
sorted_neural_data = load( path_to_sorted );
sorted_neural_data = sorted_neural_data.sorted_neural_data;
disp('Loaded');

%%

params.use_calculated_features = true;
% For dictator game manual curation, we did not have an 'unsure' category
params.include_unsure_units = true;
eisg.ei_class.do_parametric_ei_classification( sorted_neural_data, params );

%% Plot E vs I cell wfs per region

params.ei_data_filename_prefix = 'ei_classification_social_gaze';
params.default_ei_features = [1, 2];
params.default_ei_clusters = 2;
eisg.plotting.plot_e_i_wfs_by_region( sorted_neural_data, params );

%% Fetch behavioral data

path_to_trial_data = fullfile( params.processed_data_folder, params.trial_data_filename );
disp('Loading trial events data');
load( path_to_trial_data );
disp('Loaded');

%% Extract valid unit PSTH

% Remove 'everywhere' looking events
roi_ind = findnone( events.labels, 'everywhere' );
tmp_params = params;
tmp_params.mask = roi_ind;

valid_unit_psth = eisg.psth_extr.extract_valid_unit_psth( sorted_neural_data, events, tmp_params );
psth_filename = 'valid_unit_psth_social_gaze.mat';
filepath = fullfile(params.processed_data_folder, psth_filename);
disp('Saving valid unit psth');
save(filepath, 'valid_unit_psth', '-v7.3');
disp('Saved');

%% Load PSTH

psth_filename = 'valid_unit_psth_dictator_game.mat';
filepath_valid_unit_psth = fullfile(params.processed_data_folder, psth_filename);
disp('Loading valid unit PSTH');
valid_unit_psth = load(filepath_valid_unit_psth);
valid_unit_psth = valid_unit_psth.valid_unit_psth;
disp('Loaded');

%% Plot PSTH for each unit of each region and celltype

eisg.plotting.plot_psth_for_each_unit( valid_unit_psth, params );

%%




