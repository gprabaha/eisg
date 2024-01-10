%% Workspace Cleanup
clc;
clear;

%% Script Parameters
do_psth_extraction = false;
recalculate_auc = false;
smoothen_psth = true; % This is a nontrivial thing
print_stats = true;

% Validity filter for units
validity_filter = {'valid-unit', 'maybe-valid-unit'};

% For plotting
excluded_categories = {'outlier', 'ofc', 'dmpfc'};

% AUC/unit subplots
unit_auc_comparison_subplots = {'pre', 'post', 'total time'};
pre_time_range = [-0.5 0]; % in seconds
post_time_range = [0 0.5];

%% Loading Data
disp('Loading data...')
data_p = fullfile( eisg.util.project_path, 'processed_data');

% Neural data
sorted = shared_utils.io.fload( fullfile(data_p,...
  'sorted_neural_data_social_gaze.mat') );

% Behavioral data
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );

% Add task-relevant constructed ROI labels 
events = eisg.util.add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

% Celltype labels
ct_labels = shared_utils.io.fload(fullfile(data_p,...
    'celltype-labels_pfc-combined-class_p2v.mat'), 'ct_labels');
disp('Done');

%% Proprocessing Data
disp('Preprocessing data...');
[unit_spike_ts, unit_wfs, spike_labels] = eisg.util.linearize_sorted( sorted );
bfw.add_monk_labels( spike_labels );
[uuid_I, uuids] = findall( spike_labels, 'uuid',...
  find(spike_labels, validity_filter) );
match_I = bfw.find_combinations( ct_labels, uuids );
for i = 1:numel(uuid_I)
  if ( ~isempty(match_I{i}) )
    ct_label = cellstr( ct_labels, 'cell-type', match_I{i} );
    addsetcat( spike_labels, 'cell-type', ct_label, uuid_I{i} );
  end
end
replace( spike_labels, 'n', 'narrow');
replace( spike_labels, 'm', 'broad');
replace( spike_labels, 'b', 'outlier');
disp('Done');

%% Declare PSTH Extraction Parameters
min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;

%% Compute/Load PSTH and Associated Parameters
if do_psth_extraction
    disp('Computing PSTH...');
    rois = {'face', 'eyes_nf', 'whole_face', 'right_nonsocial_object_whole_face_matched'};
    evt_mask = find(events.labels, [{'m1'}, rois]);
    spk_mask = find(spike_labels, {'valid-unit', 'maybe-valid-unit'}); 
    [psth_matrix, psth_labels, t] = eisg.psth.compute_psth(...
        unit_spike_ts, spike_labels, spk_mask, ...
        evts, events.labels, evt_mask, ...
        min_t, max_t, bin_width);
    % Save variables
    disp('Saving PSTH...');
    save(fullfile(data_p, 'binned_unit_psth_social_gaze.mat')...
        , 'psth_matrix', 'psth_labels', 't'...
        );
    disp('Done');
else
    disp('Loading saved PSTH...');
    loaded_data = load(fullfile(data_p, 'binned_unit_psth_social_gaze.mat'));
    psth_matrix = loaded_data.psth_matrix;
    psth_labels = loaded_data.psth_labels;
    t = loaded_data.t;
    disp('Done');
end

%% Smoothen Data
if smoothen_psth
    disp('Smoothening PSTH...'); 
    psth_matrix = eisg.psth.smoothen_psth(psth_matrix);
else
    disp('Using raw PSTH without smoothening...');
end

%% Face vs Obj AUC
% Calculating the AUC values
roi_a = 'whole_face';
roi_b = 'right_nonsocial_object_whole_face_matched';
if recalculate_auc
    [auc_f_o, z_scored_aucs_f_o, auc_labels_f_o] = eisg.auc.calculate_roi_comparison_auc(...
        psth_matrix, psth_labels, roi_a, roi_b...
        );
    save( fullfile(data_p, 'face_obj_comparison_auc_values.mat'), 'auc_f_o', 'z_scored_aucs_f_o', 'auc_labels_f_o' );
else
    load( fullfile(data_p, 'face_obj_comparison_auc_values.mat') );
end

%% Plot and Compare AUC
% Analyze for ACC
region = 'acc';

% Ploting Face vs Obj AUC Timeseries
figure();
eisg.plot.plot_zscored_auc_timeseries_across_celltypes(...
    t, z_scored_aucs_f_o, auc_labels_f_o, region, excluded_categories);
sgtitle(['Face vs Obj AUC Timeseries for: ' region]); drawnow;

% Plotting mean AUC/unit violinplots
figure();
eisg.plot.violinplot_compare_mean_auc_per_unit(t, z_scored_aucs_f_o,...
    auc_labels_f_o, region, excluded_categories,...
    unit_auc_comparison_subplots, pre_time_range, post_time_range, print_stats);
sgtitle(['Face vs Obj AUC/unit comparison across celltypes for: ' region]); drawnow;

% Analyze for BLA
region = 'bla';

% Ploting Face vs Obj AUC Timeseries
figure();
eisg.plot.plot_zscored_auc_timeseries_across_celltypes(...
    t, z_scored_aucs_f_o, auc_labels_f_o, region, excluded_categories);
sgtitle(['Face vs Obj AUC Timeseries for: ' region]); drawnow;

% Plotting mean AUC/unit violinplots
figure();
eisg.plot.violinplot_compare_mean_auc_per_unit(t, z_scored_aucs_f_o,...
    auc_labels_f_o, region, excluded_categories,...
    unit_auc_comparison_subplots, pre_time_range, post_time_range, print_stats);
sgtitle(['Face vs Obj AUC/unit comparison across celltypes for: ' region]); drawnow;

%% EyeNF vs Face AUC
% Calculating the AUC values
roi_a = 'eyes_nf';
roi_b = 'face';
if recalculate_auc
    [auc_enf_f, z_scored_aucs_enf_f, auc_labels_enf_f] = eisg.auc.calculate_roi_comparison_auc(...
        psth_matrix, psth_labels, roi_a, roi_b...
        );
    save( fullfile(data_p, 'eyenf_face_comparison_auc_values.mat'), 'auc_enf_f', 'z_scored_aucs_enf_f', 'auc_labels_enf_f' );
else
    load( fullfile(data_p, 'eyenf_face_comparison_auc_values.mat') );
end

%% Plot and Compare AUC
% Analyze for ACC
region = 'acc';

% Ploting EyeNF vs Face AUC Timeseries
figure();
eisg.plot.plot_zscored_auc_timeseries_across_celltypes(...
    t, z_scored_aucs_enf_f, auc_labels_enf_f, region, excluded_categories);
sgtitle(['EyeNF vs Face AUC Timeseries for: ' region]); drawnow;

% Plotting mean AUC/unit violinplots
figure();
eisg.plot.violinplot_compare_mean_auc_per_unit(t, z_scored_aucs_enf_f,...
    auc_labels_enf_f, region, excluded_categories,...
    unit_auc_comparison_subplots, pre_time_range, post_time_range, print_stats);
sgtitle(['EyeNF vs Face AUC/unit comparison across celltypes for: ' region]); drawnow;

% Analyze for BLA
region = 'bla';

% Ploting EyeNF vs Face AUC Timeseries
figure();
eisg.plot.plot_zscored_auc_timeseries_across_celltypes(...
    t, z_scored_aucs_enf_f, auc_labels_enf_f, region, excluded_categories);
sgtitle(['EyeNF vs Face AUC Timeseries for: ' region]); drawnow;

% Plotting mean AUC/unit violinplots
figure();
eisg.plot.violinplot_compare_mean_auc_per_unit(t, z_scored_aucs_enf_f,...
    auc_labels_enf_f, region, excluded_categories,...
    unit_auc_comparison_subplots, pre_time_range, post_time_range, print_stats);
sgtitle(['EyeNF vs Face AUC/unit comparison across celltypes for: ' region]); drawnow;
