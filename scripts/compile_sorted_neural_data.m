
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




