%% Command Window Cleanup
clc;

%% Script Parameters
data_p             = fullfile(eisg.util.project_path, 'processed_data');
raw_behavior_root  = fullfile(data_p,            'raw_behavior');
pos_dir            = fullfile(raw_behavior_root, 'aligned_raw_samples/position');
time_dir           = fullfile(raw_behavior_root, 'aligned_raw_samples/time');
bounds_dir         = fullfile(raw_behavior_root, 'aligned_raw_samples/bounds');
fix_dir            = fullfile(raw_behavior_root, 'aligned_raw_samples/raw_eye_mmv_fixations');
meta_dir           = fullfile(raw_behavior_root, 'meta');
events_dir         = fullfile(raw_behavior_root, 'raw_events_remade');
roi_dir            = fullfile(raw_behavior_root, 'rois');

do_psth_extraction = true;

% Validity filter for units
validity_filter = {'valid-unit', 'maybe-valid-unit'};

% For analysis
excluded_categories = {'outlier', '<cell-type>', 'ofc', 'dmpfc'};


%% Loading Data
disp( 'Loading data...' );

% Neural data
if ~exist('sorted', 'var')
    disp('Loading sorted neural data...');
    sorted = shared_utils.io.fload( fullfile(data_p,...
      'sorted_neural_data_social_gaze.mat') );
else
    disp( 'Using existing sorted data in workspace' );
end

% Processed Behavioral data
if ~exist( 'events', 'var' )
    disp( 'Loading behavioral data...' );
    events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );
    % Add task-relevant constructed ROI labels 
    events = eisg.util.add_whole_face_whole_object_rois( events );
    event_start_times = bfw.event_column( events, 'start_time' );
else
    disp('Using existing behavioral events in workspace');
end

% Celltype labels
if ~exist( 'ct_labels', 'var' )
    disp( 'Loading celltype labels data...');
    ct_labels = shared_utils.io.fload(fullfile( data_p,...
        'celltype-labels_pfc-combined-class_p2v.mat'), 'ct_labels' );
else
    disp( 'Using existing ct_labels in workspace' );
end
disp( 'Done' );

%% Preprocessing Neural data

disp( 'Preprocessing data...' );
[unit_spike_ts, unit_wfs, spike_labels] = eisg.util.linearize_sorted( sorted );
bfw.add_monk_labels( spike_labels );
[uuid_I, uuids] = findall( spike_labels, 'uuid',...
  find( spike_labels, validity_filter ) );
match_I = bfw.find_combinations( ct_labels, uuids );
for i = 1:numel(uuid_I)
  if ( ~isempty( match_I{i} ) )
    ct_label = cellstr( ct_labels, 'cell-type', match_I{i} );
    addsetcat( spike_labels, 'cell-type', ct_label, uuid_I{i} );
  end
end
replace( spike_labels, 'n', 'narrow' );
replace( spike_labels, 'm', 'broad' );
replace( spike_labels, 'b', 'outlier' );
disp( 'Done' );

%% Compute event distances

pos_file_list = shared_utils.io.findmat( pos_dir );
fix_file_list = shared_utils.io.findmat( fix_dir );
roi_file_list = shared_utils.io.findmat( roi_dir );

pos_file_list(is_hidden(pos_file_list)) = [];
fix_file_list(is_hidden(fix_file_list)) = [];
roi_file_list(is_hidden(roi_file_list)) = [];
fnames_pos = shared_utils.io.filenames( pos_file_list );
fnames_roi = shared_utils.io.filenames( roi_file_list );
roi_file_list(~ismember(fnames_roi, fnames_pos)) = [];

disp('Computing fixation distances from eyes...')
event_distances = compute_event_distances_from_eyes(events, pos_file_list ...
    , fix_file_list, roi_file_list);


%% Convert to trial table

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

mask = pipe(rowmask(spike_labels) ...
    , @(m) find(spike_labels, validity_filter, m) ...
    , @(m) findnone(spike_labels, excluded_categories, m) ...
);

filtered_unit_spike_ts = unit_spike_ts(mask);
filtered_spike_labels = prune(spike_labels(mask));

spike_file = struct( ...
    'spike_times', {filtered_unit_spike_ts}, 'labels', prune(filtered_spike_labels) ...
);

%% Compute psth for all cells

sessions = trial_table.sessions;
% Celltype labels
if do_psth_extraction
    disp('Extracting PSTH')
    [spike_counts, ind_fix, ind_spikes] = compute_psth_for_gaze_tracking( ...
        spike_file, sessions, trial_table.fixation_start_ts, trial_table.fixation_stop_ts );
end

%%  convert to table

spike_count_tbl = table( ...
    categorical(filtered_spike_labels, 'uuid', ind_spikes) ...
    , categorical(filtered_spike_labels, 'region', ind_spikes) ...
    , categorical(filtered_spike_labels, 'cell-type', ind_spikes) ...
    , 'va', {'uuid', 'region', 'cell-type'} ...
);
spike_count_tbl = [ spike_count_tbl, trial_table(ind_fix, :) ];
spike_count_tbl.spike_counts = spike_counts;
%   @TODO: add contra/ipsi label for fixations
% spike_count_tbl.is_contra = spike_count_tbl.m1_hemifield_origin_delta_pos(:, 1) > 0;
spike_count_tbl.is_contra(:) = false;

%%  extract model inputs + create fit function

disp('Fitting GLM...')
[X, y, offset, base_mask] = extract_model_inputs( spike_count_tbl );
do_fit = @(ind) fit_non_step_wise_distance_model(X(ind, :), y(ind), offset(ind));

%%  run model for each cell, contra or ipsi

do_nmatch = false;

mask = base_mask;

% mask = mask & spike_count_tbl.is_contra;  % contra
% mask = mask & ~spike_count_tbl.is_contra; % ipsi

if ( do_nmatch )
  I = findeach( spike_count_tbl, 'uuid' );
  n_match_mask = n_match_contra_ipsi( spike_count_tbl.is_contra, I );
  mask = mask & n_match_mask;
end

[cell_I, mdl_tbl] = findeach( spike_count_tbl, {'uuid', 'sessions', 'region', 'cell-type'}, mask );
mdl_tbl.mdls = cell( numel(cell_I), 1 );

for i = 1:numel(cell_I)
  fprintf( '\n %d of %d', i, numel(cell_I) );
  mdl_tbl.mdls{i} = do_fit( cell_I{i} );
end

%% now iterate through each cell and create a sigcell array and label

sig_cells_m1_dist = zeros( size(mdl_tbl, 1), 1 );
sig_cells_m2_dist = zeros( size(mdl_tbl, 1), 1 );
for i=1:size(mdl_tbl, 1)
    glm = mdl_tbl.("mdls"){i};
    p_vals = glm.Coefficients.pValue;
    if p_vals(2) < 0.05
        sig_cells_m1_dist(i) = 1;
    end
    if p_vals(3) < 0.05
        sig_cells_m2_dist(i) = 1;
    end
end

mdl_tbl = addvars(mdl_tbl, sig_cells_m1_dist, sig_cells_m2_dist, 'NewVariableNames', {'sigCellsM1Dist', 'sigCellsM2Dist'});


%% Display results

% Display the fraction of sig_cells for each cell type of each region for M1 distance tracking
display_sig_cell_fractions(mdl_tbl, sig_cells_m1_dist);

plot_sig_cell_fractions_piecharts(mdl_tbl, sig_cells_m1_dist, 'M1 distance tracking');

% Perform the chi-square test for M1 distance tracking
[tbl_m1_dist, ~] = do_prop_tests(mdl_tbl, sig_cells_m1_dist);
% Display the results for M1 distance tracking
disp('Chi-Square Test for M1 Distance Tracking:');
disp(tbl_m1_dist);

% Display the fraction of sig_cells for each cell type of each region for M2 distance tracking
display_sig_cell_fractions(mdl_tbl, sig_cells_m2_dist);

plot_sig_cell_fractions_piecharts(mdl_tbl, sig_cells_m2_dist, 'M2 distance tracking');

% Perform the chi-square test for M2 distance tracking
[tbl_m2_dist, ~] = do_prop_tests(mdl_tbl, sig_cells_m2_dist);
% Display the results for M2 distance tracking
disp('Chi-Square Test for M2 Distance Tracking:');
disp(tbl_m2_dist);

%%  run model for each cell, resampling subset (contra or ipsi) with higher N

%{
num_resample = 20;

% fit for each cell
mask = base_mask;
[cell_I, mdl_tbl] = findeach( ...
  spike_count_tbl, {'uuid', 'sessions', 'region'}, mask );
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
%}


%%  Helper Functions

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


function tf = is_hidden(f)
    fnames = shared_utils.io.filenames( f );
    tf = startsWith( fnames, '.' );
end

function display_sig_cell_fractions(mdl_tbl, sig_cells)
    % Display the fraction of sig_cells for each cell type of each region
    unique_regions = unique(mdl_tbl.region);
    for r = 1:numel(unique_regions)
        disp(['Region: ', char(unique_regions(r))]);
        region_mask = mdl_tbl.region == unique_regions(r);
        unique_cell_types = unique(mdl_tbl.('cell-type')(region_mask));
        for c = 1:numel(unique_cell_types)
            cell_type_mask = mdl_tbl.('cell-type') == unique_cell_types(c) & region_mask;
            total_count = nnz(cell_type_mask);
            sig_cell_count = nnz(sig_cells(cell_type_mask));
            fraction_sig_cells = sig_cell_count / total_count;
            disp(['   Cell Type: ', char(unique_cell_types(c)), ...
                  ', Fraction of sig_cells: ', num2str(fraction_sig_cells * 100), '%', ...
                  ', Count: ', num2str(sig_cell_count), ...
                  ', Total Count: ', num2str(total_count)]);
        end
    end
end


function plot_sig_cell_fractions_piecharts(mdl_tbl, sig_cells, title_str)
    % Get unique regions
    unique_regions = unique(mdl_tbl.region);
    
    % Loop through each region
    for r = 1:numel(unique_regions)
        % Extract region-specific data
        region_name = char(unique_regions(r));
        region_mask = mdl_tbl.region == unique_regions(r);
        
        % Get unique cell types within the region
        unique_cell_types = unique(mdl_tbl.('cell-type')(region_mask));
        
        % Initialize figure and subplots for the current region
        figure;
        sgtitle(['Sig Cell Fractions for Region: ', region_name, ' - ', title_str]);
        subplot_cols = 2;  % Set the number of subplot columns
        subplot_rows = ceil(numel(unique_cell_types) / subplot_cols);
        
        % Loop through each cell type
        for c = 1:numel(unique_cell_types)
            % Extract cell type-specific data
            cell_type_name = char(unique_cell_types(c));
            cell_type_mask = mdl_tbl.('cell-type') == unique_cell_types(c) & region_mask;
            
            % Calculate sig cell fraction
            total_count = nnz(cell_type_mask);
            sig_cell_count = nnz(sig_cells(cell_type_mask));
            fraction_sig_cells = sig_cell_count / total_count;
            
            % Convert sig_cell_count and total_count to numerical arrays
            count_array = [sig_cell_count, total_count - sig_cell_count];
            
            % Create subplot for the current cell type
            subplot(subplot_rows, subplot_cols, c);
            
            % Define explode vector to highlight sig_cells
            explode = [0.1, 0];
            % Plot pie chart with explode and labels
            h = pie(count_array, explode, {'Sig Cells', 'Non-Sig Cells'});
            % Set color of sig_cells slice to red
            h(1).FaceColor = 'r';
            % Add title
            title([cell_type_name, ' - ', title_str]);
            
            % Display total number of sig_cells and percentage below the pie chart
            total_str = sprintf('Sig Cell Fraction: %.2f%% (%d/%d total units)', fraction_sig_cells * 100, sig_cell_count, total_count);
            text(0, -1.2, total_str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12);
        end
    end
end


function [tbl, reg_labels, tbls] = do_prop_tests(mdl_tbl, sig_cells)
    % Extract relevant columns from mdl_tbl
    regions = mdl_tbl.region;
    cell_types = mdl_tbl.('cell-type');

    % Initialize variables to store results
    tbls = struct('h', {}, 'p', {}, 'chi2stat', {});

    % Loop through each unique region
    unique_regions = unique(regions);
    for r = 1:numel(unique_regions)
        % Find rows corresponding to the current region
        region_mask = regions == unique_regions(r);
        
        % Get unique cell types within the current region
        unique_cell_types = unique(cell_types(region_mask));
        
        % Initialize variables to store counts for each cell type
        sig_cell_counts = zeros(1, numel(unique_cell_types));
        total_counts = zeros(1, numel(unique_cell_types));
        
        % Loop through each unique cell type within the current region
        for c = 1:numel(unique_cell_types)
            % Find rows corresponding to the current cell type
            cell_type_mask = cell_types == unique_cell_types(c);
            
            % Calculate counts for the current cell type
            sig_cell_counts(c) = nnz(sig_cells(region_mask & cell_type_mask));
            total_counts(c) = nnz(region_mask & cell_type_mask);
        end
        
        % Perform the proportion test for pairwise comparisons of cell types
        [h, p, chi2stat] = eisg.util.prop_test(sig_cell_counts, total_counts, false);
        
        % Store the results
        tbls(r) = struct('h', h, 'p', p, 'chi2stat', chi2stat);
    end

    % Convert results to table format
    reg_labels = unique_regions'; % Replace spaces with underscores
    tbl = array2table([tbls.p]', 'RowNames', cellstr(reg_labels), 'VariableNames', {'Chi2-pValue'});
end

%{
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
%}