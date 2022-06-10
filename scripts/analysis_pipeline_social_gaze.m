
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
  'valley_ind', ...
  'peak_ind', ...
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
params.include_unsure_units = true;

%% Loading sorted neural data

path_to_sorted = fullfile( params.processed_data_folder, params.sorted_data_filename );
disp('Loading sorted neural data');
sorted_neural_data = load( path_to_sorted );
sorted_neural_data = sorted_neural_data.sorted_neural_data;
disp('Loaded');

%%

valid_unit_features = eisg.util.extract_valid_unit_features_from_sorted(sorted_neural_data, params);

%% Classify celltypes based on peak-to-valley

regions = unique({sorted_neural_data.region});
regions = regions([1,2,3]);

n_cluster_range = [2 7];

% 'b' stands for 'broad, 'm' for 'medium', and 'n' for 'narrow' spiking
% units
celltype_labels = eisg.celltype_class.do_xiang_multi_cluster_kmeans_by_region( valid_unit_features, regions, n_cluster_range );

label_filename = 'celltype_labels_classified_by_region.mat';
label_path = fullfile( params.processed_data_folder, label_filename );
save(label_path, 'celltype_labels');

%% Fetch behavioral data

params.trial_data_filename = 'trial_events_social_gaze.mat';
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

psth_filename = 'valid_unit_psth_social_gaze.mat';
filepath_valid_unit_psth = fullfile(params.processed_data_folder, psth_filename);
disp('Loading valid unit PSTH');
valid_unit_psth = load(filepath_valid_unit_psth);
valid_unit_psth = valid_unit_psth.valid_unit_psth;
disp('Loaded');

%% Plot PSTH for each unit, separated by region and celltype


