function params = make_analysis_params_struct()
 
params = struct();
 
% Common, and PSTH extraction
params.processed_data_folder = 'processed_data';
params.sorted_data_filename = 'sorted_neural_data_social_gaze.mat';
params.trial_data_filename = [];
params.include_unsure_units = false;
params.sorted_data_filename = [];
params.trial_data_filename = [];
params.ei_label_filename = [];

% WF features
params.feature_list = {
    'peak_to_valley', ...
    'halfwidth', ...
    'peak_trough_ratio', ...
    'recovery_slope', ...
    'repolarization_slope'
};

% Celltype class
params.use_calculated_features = true;

% ROI-s
params.face_roi = {'face', 'eyes_nf'};
params.eye_roi = {'eyes_nf'};
params.ne_face_roi = {'face'};
params.obj_roi = {'left_nonsocial_object', 'right_nonsocial_object'};

% Analysis preiods
params.pre_fix_time_win = [-0.4 -0.05];
params.post_fix_time_win = [0.05 0.4];

% PSTH plotting and analysis
params.plot_folder = 'plots';
 
params.region = 'bla';
params.tasktype = 'free_viewing';
params.population_psth_method = 'mean';

 
% Progress bar
params.num_ticks_in_progress_bar = 75;
 
end
