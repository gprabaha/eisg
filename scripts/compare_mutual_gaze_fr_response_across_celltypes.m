%% Command Window Cleanup
clc;

%% Script Parameters
script_params = struct();
script_params.data_p                    = fullfile( eisg.util.project_path, 'processed_data' );
% Unit selection
script_params.validity_filter           = {'valid-unit', 'maybe-valid-unit'};
% Comparison method
script_params.method                    = 'ranksum'; % available: ranksum, ttest2
% Time window for mean/median calculation
script_params.analysis_time_window      = [-0.2 0.2]; % [-0.2 0.2] is pretty much the only time window that works
script_params.excluded_categories       = {'outlier', 'ofc', 'dmpfc'};
% PSTH related
script_params.do_psth_extraction        = false;
script_params.smoothen_psth             = true; % smoothening gives better results for [-0.2 to 0.2]
script_params.min_t                     = -0.5;
script_params.max_t                     = 0.5;
script_params.bin_width                 = 0.01;
% Specific to mutual gaze
script_params.roi                       = 'eyes_nf'; % analyze only looks into eyes_nf
script_params.mutual_associated_labels  = {'m1', 'm2', 'm1_initiated', 'm2_initiated', 'mutual'};

disp( 'Parameter values used for this script run are:' );
disp( script_params );

% !!
% ranksum, with unsmoothened psth in the [-0.25 to 0.25]s time window shows
% interesting results!
% !!

%% Loading Data
data_p = script_params.data_p;
disp( 'Loading data...' );

% Neural data
if ~exist('sorted', 'var')
    disp('Loading sorted neural data...');
    sorted = shared_utils.io.fload( fullfile(data_p,...
      'sorted_neural_data_social_gaze.mat') );
else
    disp( 'Using existing sorted data in workspace' );
end

% Behavioral data
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

%% Proprocessing Data
disp( 'Preprocessing data...' );
[unit_spike_ts, unit_wfs, spike_labels] = eisg.util.linearize_sorted( sorted );
bfw.add_monk_labels( spike_labels );
[uuid_I, uuids] = findall( spike_labels, 'uuid',...
  find( spike_labels, script_params.validity_filter ) );
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

%% Compute/Load PSTH
do_psth_extraction = script_params.do_psth_extraction;
if do_psth_extraction
    disp( 'Computing PSTH...' );
    rois = script_params.roi;
    mutual_associated_labels = script_params.mutual_associated_labels;
    validity_filter = script_params.validity_filter;
    evt_mask = find(events.labels, [mutual_associated_labels, rois]);
    spk_mask = find(spike_labels, validity_filter); 
    [psth_matrix, psth_labels, t] = eisg.psth.compute_psth(...
        unit_spike_ts, spike_labels, spk_mask, ...
        event_start_times, events.labels, evt_mask, ...
        script_params.min_t, script_params.max_t, script_params.bin_width);
    % Save variables
    disp( 'Saving PSTH...' );
    save( fullfile( data_p, 'binned_eye_nf_psth_mutual_gaze.mat' )...
        , 'psth_matrix', 'psth_labels', 't'...
        );
    disp( 'Done' );
elseif exist( 'psth_matrix', 'var' )
    disp( 'Using existing PSTH matrix in workspace' );
else
    disp( 'Loading saved PSTH...' );
    loaded_data = load( fullfile( data_p, 'binned_eye_nf_psth_mutual_gaze.mat' ) );
    psth_matrix = loaded_data.psth_matrix;
    psth_labels = loaded_data.psth_labels;
    t = loaded_data.t;
    disp( 'Done' );
end

%% Smoothen Data
smoothen_psth = script_params.smoothen_psth;
if smoothen_psth
    disp( 'Smoothening PSTH... ' ); 
    spike_ct_mat = eisg.psth.smoothen_psth( psth_matrix );
    disp( 'Done' );
else
    disp( 'Using raw PSTH without smoothening' );
    spike_ct_mat = psth_matrix;
end

%% Compare FR Response
excluded_categories = script_params.excluded_categories;
mask = findnone( psth_labels, excluded_categories );
unit_labels = retaineach( psth_labels, {'uuid'}, mask ); % Select 1 row of labels for each uuid
[I, C] = findall( psth_labels, {'uuid'}, mask );
sig_units_m1_mutual_vs_excl             = nan( size(unit_labels, 1), 1 );
sig_units_m2_mutual_vs_excl             = nan( size(unit_labels, 1), 1 );
sig_units_mutual_vs_m1_excl             = nan( size(unit_labels, 1), 1 );
sig_units_m1_init_vs_m2_init_mutual     = nan( size(unit_labels, 1), 1 );
for i = 1:numel(I)
    inds = I{i};
    roi = script_params.roi;
    m1_excl_inds = find( psth_labels, {'m1', roi}, inds );
    m1_init_inds = find( psth_labels, {'mutual', 'm1_initiated', roi}, inds );
    m2_excl_inds = find( psth_labels, {'m2', roi}, inds );
    m2_init_inds = find( psth_labels, {'mutual', 'm2_initiated', roi}, inds );
    method = script_params.method;
    analysis_time_window = script_params.analysis_time_window;
    sig_units_m1_mutual_vs_excl(i) = eisg.mean_fr.stat_compare_mean_fr(...
        spike_ct_mat, t, m1_excl_inds, m1_init_inds, method, analysis_time_window );
    sig_units_m2_mutual_vs_excl(i) = eisg.mean_fr.stat_compare_mean_fr(...
        spike_ct_mat, t, m2_excl_inds, m2_init_inds, method, analysis_time_window );
    sig_units_mutual_vs_m1_excl(i) = eisg.mean_fr.stat_compare_mean_fr(...
        spike_ct_mat, t, m1_excl_inds, [m1_init_inds; m2_init_inds], method, analysis_time_window );
    sig_units_m1_init_vs_m2_init_mutual(i) = eisg.mean_fr.stat_compare_mean_fr(...
        spike_ct_mat, t, m1_init_inds, m2_init_inds, method, analysis_time_window );
end

%% Celltype proportion comparison
% M1 initated mutual vs exclusive
display_celltype_prop_stats( sig_units_m1_mutual_vs_excl, unit_labels,...
    '** M1 initated mutual vs exclusive events **' );

% M2 initated mutual vs exclusive
display_celltype_prop_stats( sig_units_m2_mutual_vs_excl, unit_labels,...
    '** M2 initated mutual vs exclusive events **' );

% All mutual vs M1 exclusive events
display_celltype_prop_stats( sig_units_mutual_vs_m1_excl, unit_labels,...
    '** All mutual vs M1 exclusive events **' );

% M1 initiated vs M2 initiated mutual events
display_celltype_prop_stats( sig_units_m1_init_vs_m2_init_mutual, unit_labels,...
    '** M1 initiated vs M2 initiated mutual events **' );

%% Helper Functions
function display_celltype_prop_stats(sig_unit_vec, unit_labels, display_string)
    sig_units = sig_unit_vec;
    %{
    sig_units(isnan( sig_units )) = 0;
    use_labels = unit_labels;
    %}
    % %{
    sig_units(isnan( sig_units )) = [];
    use_labels = unit_labels(find( ~isnan( sig_units ) ));
    % %}
    mask = rowmask( use_labels );
    [~, unit_inds_in_comb, celltype_region_comb] = keepeach( use_labels', {'cell-type', 'region'} );
    cell_type_count = cellfun( @(x) numel( sig_units(x) ), unit_inds_in_comb );
    sig_cell_count = cellfun( @(x) nnz( sig_units(x) ), unit_inds_in_comb );
    sig_cell_frac = cellfun( @(x) pnz( sig_units(x) ), unit_inds_in_comb );
    prop_table = table( cell_type_count, sig_cell_count, sig_cell_frac, 'RowNames', fcat.strjoin(celltype_region_comb, '_') );
    tbl = do_prop_tests( sig_units, use_labels', mask );

    disp( display_string );
    disp( repmat( '-', 1, numel( display_string ) ) );
    disp( prop_table );
    disp( 'Chi-square test comparing celltype proportions:' );
    disp( tbl );
end

function [tbl, reg_labels, tbls] = do_prop_tests(meets_criterion, labels, mask)
    assert_ispair(meets_criterion, labels);

    [reg_labels, reg_I] = keepeach(labels', 'region', mask);
    n_I = cellfun(@(x) find(labels, 'narrow', x), reg_I, 'un', 0);
    m_I = cellfun(@(x) find(labels, 'broad', x), reg_I, 'un', 0);

    tbls = do_prop_test(meets_criterion, n_I, m_I);
    chi2ps = [tbls.p]';
    tbl = array2table(chi2ps, 'RowNames', fcat.strjoin(reg_labels(:, 'region')'));
end

function tbls = do_prop_test(meets_criterion, inds_a, inds_b)
    xa = cellfun(@(x) sum(meets_criterion(x)), inds_a);
    xb = cellfun(@(x) sum(meets_criterion(x)), inds_b);
    na = cellfun(@numel, inds_a);
    nb = cellfun(@numel, inds_b);

    tbls = struct('h', {}, 'p', {}, 'chi2stat', {});
    for i = 1:numel(xa)
        [h, p, chi2stat] = eisg.util.prop_test([xa(i), xb(i)], [na(i), nb(i)], false);
        tbls(i) = struct('h', h, 'p', p, 'chi2stat', chi2stat);
    end
end
