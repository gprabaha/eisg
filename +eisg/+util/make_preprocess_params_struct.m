function params = make_preprocess_params_struct()

params = struct();

params.processed_data_folder = 'processed_data';
params.sorted_data_filename = 'sorted_neural_data_social_gaze.mat';
params.include_unsure_units = true;
params.path_to_sorter_output = [];
params.sorted_neural_data_filename = [];
params.path_to_validity = [];
params.integrate_validity = true;
params.calculate_features = true;
params.feature_list = [];
params.n_extra_pts_per_gap_for_interpolation = 2;
params.num_ticks_in_progress_bar = 75;

end