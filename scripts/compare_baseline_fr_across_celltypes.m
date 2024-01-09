%% Workspace Cleanup
clc;
clear;

%% Script Parameters
regenerate_psth_for_all_roi = false;
recalculate_baseline_stats = false;

%% Loading Data
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

%% Data Preprocessing
disp('Preprocessing data...');

% Extract 2D fcat matrix and labels from spike struct and unit WFs
[unit_spike_ts, unit_wfs, spike_labels] = eisg.util.linearize_sorted(sorted);
bfw.add_monk_labels(spike_labels); % Add labels of participant monkeys

% Update cell-type Labels
disp('Updating cell-type labels...');
replace(ct_labels, 'n', 'narrow');
replace(ct_labels, 'm', 'broad');
replace(ct_labels, 'b', 'outlier');

%% Parameters to Extract Mean fr from PSTH
disp('Declaring parameters for PSTH extraction...');
min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;

% Validity filter for units
validity_filter = {'valid-unit', 'maybe-valid-unit'};
all_fix = {'everywhere', 'right_nonsocial_object_whole_face_matched', 'whole_face'};

%% PSTH Extraction Across All ROIs
% Really heavy to generate on the RAM and takes ample time
if regenerate_psth_for_all_roi
    disp('Extracting PSTH...');
    evt_mask = find(events.labels, [{'m1'}, all_fix]);
    spk_mask = find(spike_labels, validity_filter);
    [psth_matrix, psth_labels, t] = eisg.psth.compute_psth(...
        unit_spike_ts, spike_labels, spk_mask, ...
        evts, events.labels, evt_mask, ...
        min_t, max_t, bin_width);
    disp('Saving output...')
    save(fullfile(data_p, 'binned_unit_psth_social_gaze_all_fix.mat'),...
        'psth_matrix', 'psth_labels', 't'...
        , '-v7.3');
else
    disp('Loading previously saved PSTH...');
    load(fullfile(data_p, 'binned_unit_psth_social_gaze_all_fix.mat'));
end

%% Declare Baseline Parameters
baseline_time_snippet = [-0.4, -0.05];
sprintf('Baseline time snippet is from %0.2f s to %0.2f s\n', baseline_time_snippet(1), baseline_time_snippet(2));

%% Calculate Baseline Stats
if recalculate_baseline_stats
    disp('Recalculating unit baseline stats...');
    
    % Baseline mean miring across all fixations for each unit
    fixation_query = all_fix;
    fixation_mask = find(psth_labels, fixation_query);
    disp('Calculating stats for all fixations...');
    [unit_I, uuids] = findall(psth_labels, 'uuid', fixation_mask);
    [mean_fr_all_fix, mean_fr_labs_all_fix] = ...
        eisg.mean_fr.unit_baseline_fr_across_fix(...
        psth_matrix, psth_labels, t, unit_I, baseline_time_snippet...
        );
    
    % Baseline mean firing across non-roi fixations for each unit
    fixation_query = {'everywhere'};
    fixation_mask = find(psth_labels, fixation_query);
    disp('Calculating stats for non-ROI fixations...');
    [unit_I, uuids] = findall(psth_labels, 'uuid', fixation_mask);
    [mean_fr_non_roi_fix, mean_fr_labs_non_roi_fix] = ...
        eisg.mean_fr.unit_baseline_fr_across_fix(...
        psth_matrix, psth_labels, t, unit_I, baseline_time_snippet...
        );
    
    % Save output
    disp('Saving output...')
    save(fullfile(data_p, 'mean_fr_for_each_unit.mat')...
        , 'mean_fr_all_fix', 'mean_fr_labs_all_fix'...
        , 'mean_fr_non_roi_fix', 'mean_fr_labs_non_roi_fix'...
        , 'baseline_time_snippet'...
        );
else
    disp('Loading previously saved base stats...');
    load(fullfile(data_p, 'mean_fr_for_each_unit.mat'));
end


%% Plot violinplots comparing mean fr (all fix) of celltypes across regions

% Fig saving parameters
save_p = fullfile(eisg.util.project_path, 'data/plots/celltype_mean_fr', dsp3.datedir);
do_save = false;

% Add cell-type labels to a temporary spike_labels variable
[uuid_I, uuid_C] = findall(ct_labels, {'uuid', 'cell-type'});
temp_fr_labels = mean_fr_labs_all_fix;
match_inds = bfw.find_combinations(temp_fr_labels, uuid_C(1, :));
addcat(temp_fr_labels, 'cell-type');
for i = 1:numel(match_inds)
    setcat(temp_fr_labels, 'cell-type', uuid_C{2, i}, match_inds{i});
end

mask = pipe(rowmask(temp_fr_labels) ...
    , @(m) find(temp_fr_labels, validity_filter, m) ...
    , @(m) findnone(temp_fr_labels, '<cell-type>', m) ...
);


pl = plotlabeled.make_common();
axs = pl.violinplot(mean_fr_all_fix(mask, 1), prune(temp_fr_labels(mask)), 'cell-type', 'region');

if (do_save)
    dsp3.req_savefig(gcf, save_p, prune(temp_fr_labels), 'region', ['mean_fr_all_fix_across_regional_celltypes_classified_by' is_pca_string(use_pca) get_feature_string_from_labels(feature_list)]);
end

% Check if the mean FRs are statistically distinct across cell types
firing_rate_subset = mean_fr_all_fix(mask);
temp_fr_labels = prune(temp_fr_labels(mask));
regions = temp_fr_labels('region')';

disp('T-Test | All Fixations');
for reg = regions
    broad_unit_inds = find(temp_fr_labels, {char(reg), 'broad'});
    narrow_unit_inds = find(temp_fr_labels, {char(reg), 'narrow'});
    [~, p] = ttest2(firing_rate_subset(broad_unit_inds), firing_rate_subset(narrow_unit_inds));
    fprintf('%s | Mean FR | Broad = %0.3f | Narrow = %0.3f | p-val = %0.3f\n' ...
        , char(reg) ...
        , nanmean(firing_rate_subset(broad_unit_inds)) ...
        , nanmean(firing_rate_subset(narrow_unit_inds)) ...
        , p);
end

disp('Ranksum | All Fixations');
for reg = regions
    broad_unit_inds = find(temp_fr_labels, {char(reg), 'broad'});
    narrow_unit_inds = find(temp_fr_labels, {char(reg), 'narrow'});
    [p, ~] = ranksum(firing_rate_subset(broad_unit_inds), firing_rate_subset(narrow_unit_inds));
    fprintf('%s | Median FR | Broad = %0.3f | Narrow = %0.3f | ranksum p-val = %0.3f\n' ...
        , char(reg) ...
        , nanmedian(firing_rate_subset(broad_unit_inds)) ...
        , nanmedian(firing_rate_subset(narrow_unit_inds)) ...
        , p);
end


%% Plot violinplots comparing mean fr (non roi fix) of celltypes across regions

% Fig saving parameters
save_p = fullfile(eisg.util.project_path, 'data/plots/celltype_mean_fr', dsp3.datedir);
do_save = false;

% Add cell-type labels to a temporary spike_labels variable
[uuid_I, uuid_C] = findall(ct_labels, {'uuid', 'cell-type'});
temp_fr_labels = mean_fr_labs_non_roi_fix;
match_inds = bfw.find_combinations(temp_fr_labels, uuid_C(1, :));
addcat(temp_fr_labels, 'cell-type');
for i = 1:numel(match_inds)
    setcat(temp_fr_labels, 'cell-type', uuid_C{2, i}, match_inds{i});
end

mask = pipe(rowmask(temp_fr_labels) ...
    , @(m) find(temp_fr_labels, validity_filter, m) ...
    , @(m) findnone(temp_fr_labels, '<cell-type>', m) ...
);

figure();
pl = plotlabeled.make_common();
axs = pl.violinplot(mean_fr_non_roi_fix(mask, 1), prune(temp_fr_labels(mask)), {'cell-type'}, 'region');

if (do_save)
    dsp3.req_savefig(gcf, save_p, prune(temp_fr_labels), 'region',...
        ['mean_fr_non_roi_fix_across_regional_celltypes_classified_by'...
        is_pca_string(use_pca) get_feature_string_from_labels(feature_list)]);
end

% Comparing mean FR across celltypes again but now for non-roi fixations
firing_rate_subset = mean_fr_non_roi_fix(mask);
temp_fr_labels = prune(temp_fr_labels(mask));
regions = temp_fr_labels('region')';

disp('T-Test | Non-ROI Fixations');
for reg = regions
    broad_unit_inds = find(temp_fr_labels, {char(reg), 'broad'});
    narrow_unit_inds = find(temp_fr_labels, {char(reg), 'narrow'});
    [~, p] = ttest2(firing_rate_subset(broad_unit_inds), firing_rate_subset(narrow_unit_inds));
    fprintf('%s | Mean FR | Broad = %0.3f | Narrow = %0.3f | ttest2 p-val = %0.3f\n' ...
        , char(reg) ...
        , nanmean(firing_rate_subset(broad_unit_inds)) ...
        , nanmean(firing_rate_subset(narrow_unit_inds)) ...
        , p);
end

disp('Ranksum | Non-ROI Fixations');
for reg = regions
    broad_unit_inds = find(temp_fr_labels, {char(reg), 'broad'});
    narrow_unit_inds = find(temp_fr_labels, {char(reg), 'narrow'});
    [p, ~] = ranksum(firing_rate_subset(broad_unit_inds), firing_rate_subset(narrow_unit_inds));
    fprintf('%s | Median FR | Broad = %0.3f | Narrow = %0.3f | ranksum p-val = %0.3f\n' ...
        , char(reg) ...
        , nanmedian(firing_rate_subset(broad_unit_inds)) ...
        , nanmedian(firing_rate_subset(narrow_unit_inds)) ...
        , p);
end


%% Old method of getting baseline stats

%% PSTH Extraction for All Fixations
% disp('Extracting PSTH of each unit corresponting to various gaze events...');
% rois_all_fix = { 'everywhere', 'right_nonsocial_object_whole_face_matched', 'whole_face' };
% evt_mask = find( events.labels, [{'m1'}, rois_all_fix] );
% spk_mask = rowmask( spike_labels );
% [mean_fr_all_fix, mean_fr_labs_all_fix] = baseline_psth_stats( ...
%   evts, events.labels, evt_mask, unit_spike_ts, spike_labels, spk_mask ...
%   , 'min_t', min_t ...
%   , 'max_t', max_t ...
%   , 'bin_width', bin_width ...
% );

% non_roi_fix = { 'everywhere' };
% evt_mask = find( events.labels, [{'m1'}, non_roi_fix] );
% 
% [mean_fr_non_roi_fix, mean_fr_labs_non_roi_fix] = baseline_psth_stats( ...
%   evts, events.labels, evt_mask, unit_spike_ts, spike_labels, spk_mask ...
%   , 'min_t', min_t ...
%   , 'max_t', max_t ...
%   , 'bin_width', bin_width ...
% );

