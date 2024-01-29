%%  load in data

% data_dir = '/Volumes/Brains/BRAINS Recording to Stim';
data_dir = '/Volumes/external3/data/changlab/siqi/distance_model/Data/reformatted';

trial_table = shared_utils.io.fload( fullfile(data_dir, 'trial_table.mat') );
spike_file = shared_utils.io.fload( fullfile(data_dir, 'relabeled_cells.mat') );

bfw.add_monk_labels( spike_file.labels );
spike_labels = fcat.totable( spike_file.labels );

%%

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

%%

intermediate_root = '/Volumes/ExtSSD/changlab/behavioral_data/social_gaze/';
pos_file_list = shared_utils.io.findmat( fullfile(intermediate_root, 'aligned_raw_samples/position'));
fix_file_list = shared_utils.io.findmat( fullfile(intermediate_root, 'aligned_raw_samples/raw_eye_mmv_fixations'));
roi_file_list = shared_utils.io.findmat( fullfile(intermediate_root, 'rois'));

pos_file_list(is_hidden(pos_file_list)) = [];
fix_file_list(is_hidden(fix_file_list)) = [];
roi_file_list(is_hidden(roi_file_list)) = [];
fnames_pos = shared_utils.io.filenames( pos_file_list );
fnames_roi = shared_utils.io.filenames( roi_file_list );
roi_file_list(~ismember(fnames_roi, fnames_pos)) = [];

event_distances = compute_event_distances(events, pos_file_list, fix_file_list, roi_file_list);

%% Preprocessing Data
disp('Preprocessing loaded data...')
[unit_spike_ts, unit_wfs, spike_labels] = eisg.util.linearize_sorted(sorted);
bfw.add_monk_labels(spike_labels);
[uuid_I, uuids] = findall(spike_labels, 'uuid', find(spike_labels, validity_filter));
match_I = bfw.find_combinations(ct_labels, uuids);
for i = 1:numel(uuid_I)
    if (~isempty(match_I{i}))
        ct_label = cellstr(ct_labels, 'cell-type', match_I{i});
        addsetcat(spike_labels, 'cell-type', ct_label, uuid_I{i});
    end
end
replace(spike_labels, 'n', 'narrow');
replace(spike_labels, 'm', 'broad');
replace(spike_labels, 'b', 'outlier');

%%  convert to trial table

trial_table = table( ...
    bfw.event_column(events, 'start_time') ...
    , bfw.event_column(events, 'stop_time') ...
    , categorical(events.labels, 'session') ...
    , event_distances.m1_dist_to_m2s_eyes ...
    , event_distances.m2_dist_to_m1s_eyes ...
    , event_distances.m2_fix_props >= 0.9 ...
    , 'va', {'fixation_start_ts', 'fixation_stop_ts', 'sessions' ...
    , 'distances_to_m2s_eyes', 'm2_dist_to_m1s_eyes', 'm2_dist_is_valid'} ...
);

trial_table.m1_fix_durations = trial_table.fixation_stop_ts - trial_table.fixation_start_ts;

spike_file = struct( ...
    'spike_times', {unit_spike_ts}, 'labels', spike_labels ...
);

%%  compute psth for all cells

sessions = trial_table.sessions;

[spike_counts, ind_fix, ind_spikes] = compute_psth( ...
  spike_file, sessions, trial_table.fixation_start_ts, trial_table.fixation_stop_ts );

%%  convert to table

spike_count_tbl = table( ...
    categorical(spike_labels, 'uuid', ind_spikes) ...
    , categorical(spike_labels, 'region', ind_spikes) ...
    , 'va', {'uuid', 'region'} ...
);
spike_count_tbl = [ spike_count_tbl, trial_table(ind_fix, :) ];
spike_count_tbl.spike_counts = spike_counts;
%   @TODO: add contra/ipsi label for fixations
% spike_count_tbl.is_contra = spike_count_tbl.m1_hemifield_origin_delta_pos(:, 1) > 0;
spike_count_tbl.is_contra(:) = false;

%%  extract model inputs + create fit function

[X, y, offset, base_mask] = extract_model_inputs( spike_count_tbl );
do_fit = @(ind) fit_non_step_wise_distance_model(X(ind, :), y(ind), offset(ind));

%%  run model for each cell, contra or ipsi

do_nmatch = false;

mask = base_mask;

% mask = mask & spike_count_tbl.is_contra;  % contra
% mask = mask & ~spike_count_tbl.is_contra; % ipsi

if ( do_nmatch )
  I = findeach( spike_count_tbl, 'unit_uuid' );
  n_match_mask = n_match_contra_ipsi( spike_count_tbl.is_contra, I );
  mask = mask & n_match_mask;
end

[cell_I, mdl_tbl] = findeach( spike_count_tbl, {'uuid', 'sessions', 'region'}, mask );
mdl_tbl.mdls = cell( numel(cell_I), 1 );

for i = 1:numel(cell_I)
  fprintf( '\n %d of %d', i, numel(cell_I) );
  mdl_tbl.mdls{i} = do_fit( cell_I{i} );
end

%%  run model for each cell, resampling subset (contra or ipsi) with higher N

num_resample = 20;

% fit for each cell
mask = base_mask;
[cell_I, mdl_tbl] = findeach( ...
  spike_count_tbl, {'unit_uuid', 'sessions', 'region'}, mask );
mdl_tbl.mdl_contras = cell( numel(cell_I), 1 );
mdl_tbl.mdl_ipsis = cell( numel(cell_I), 1 );

for i = 1:numel(cell_I)
  fprintf( '\n %d of %d', i, numel(cell_I) );
  
  ci = cell_I{i};
  
  ci_contra = ci(spike_count_tbl.is_contra(ci));
  ci_ipsi = ci(~spike_count_tbl.is_contra(ci));
  
  nc = numel( ci_contra );
  ni = numel( ci_ipsi );
  
  if ( ni >= nc )
    mdl_ipsis = cell( num_resample, 1 );
    mdl_contras = { do_fit(ci_contra) };
    
    for j = 1:num_resample
      resamp = randsample( ni, nc, true );
      ind = ci_ipsi(resamp);
      mdl_ipsis{j} = do_fit( ind );
    end
  else
    mdl_contras = cell( num_resample, 1 );
    mdl_ipsis = { do_fit(ci_ipsi) };
    
    for j = 1:num_resample
      resamp = randsample( nc, ni, true );
      ind = ci_contra(resamp);
      mdl_contras{j} = do_fit( ind );
    end
  end
  
  mdl_tbl.mdl_contras{i} = mdl_contras;
  mdl_tbl.mdl_ipsis{i} = mdl_ipsis;
end

%%  

function [X, y, offset, mask] = extract_model_inputs(spike_count_tbl)

self_distance = spike_count_tbl.distances_to_m2s_eyes;
self_distance(self_distance > 20) = nan;

other_distance = spike_count_tbl.m2_dist_to_m1s_eyes;
other_distance(~spike_count_tbl.m2_dist_is_valid) = nan;

% fit self + other distance
X = [self_distance, other_distance];

% only fit self distance
% X = [self_distance];

% fit to spike counts
y = spike_count_tbl.spike_counts;

% offset by fixation duration in s
offset = spike_count_tbl.m1_fix_durations / 1e3;

% keep only valid rows
mask = ~any( isnan(X), 2 ) & ~isnan( y ) & ~isnan( offset );

end

function mask = n_match_contra_ipsi(is_contra, I)

mask = false( size(is_contra) );

for i = 1:numel(I)
  ic = I{i}(is_contra(I{i}));
  ii = I{i}(~is_contra(I{i}));
  
  if ( numel(ii) > numel(ic) )
    % match ipsi to contra, keep all contra
    keep_ni = ii(randperm(numel(ii), numel(ic)));
    keep_nc = ic;
  elseif ( numel(ic) > numel(ii) )
    % match contra to ipsi, keep all ipsi
    keep_nc = ic(randperm(numel(ic), numel(ii)));
    keep_ni = ii;
  else
    keep_ni = ii;
    keep_nc = ic;
  end
  
  mask(keep_ni) = true;
  mask(keep_nc) = true;
end

end

function tf = is_hidden(f)

fnames = shared_utils.io.filenames( f );
tf = startsWith( fnames, '.' );

end