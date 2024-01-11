%% Workspace Cleanup
clc;
clear;

%% Script Parameters
%  PSTH parameters
do_psth_extraction = false;
smoothen_psth = false;
min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;

% Validity filter for units
validity_filter = {'valid-unit', 'maybe-valid-unit'};

% Mean calculation parameters
method = 'ranksum'; % Method to compare FR vectors
% method = 'ttest2'; % Method to compare FR vectors
% analysis_time_window = [0 0.5];
% analysis_time_window = [-0.5 0];
analysis_time_window = [-0.25 0.25];
% analysis_time_window = [-0.5 0.5];

% !!
% ranksum, with unsmoothened psth in the [-0.25 to 0.25]s time window shows
% interesting results!
% !!

% Data discarded
excluded_categories = {'outlier', 'ofc', 'dmpfc'};

%% Loading Data
data_p = fullfile( eisg.util.project_path, 'processed_data');
fprintf('Data folder path is: %s\n', data_p);
disp('Loading data...');

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


%% Compute/Load PSTH
if do_psth_extraction
    disp('Computing PSTH...');
    rois = {'eyes_nf'};
    mutual_associated_labels = {'m1', 'm2', 'm1_initiated', 'm2_initiated', 'mutual'};
    evt_mask = find(events.labels, [mutual_associated_labels, rois]);
    spk_mask = find(spike_labels, {'valid-unit', 'maybe-valid-unit'}); 
    [psth_matrix, psth_labels, t] = eisg.psth.compute_psth(...
        unit_spike_ts, spike_labels, spk_mask, ...
        evts, events.labels, evt_mask, ...
        min_t, max_t, bin_width);
    % Save variables
    disp('Saving PSTH...');
    save(fullfile(data_p, 'binned_eye_nf_psth_mutual_gaze.mat')...
        , 'psth_matrix', 'psth_labels', 't'...
        );
    disp('Done');
else
    disp('Loading saved PSTH...');
    loaded_data = load(fullfile(data_p, 'binned_eye_nf_psth_mutual_gaze.mat'));
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

%% Compare FR Response

mask = findnone( psth_labels, excluded_categories );
unit_labels = retaineach( psth_labels, {'uuid'}, mask ); % Select 1 row of labels for each uuid
m1_mutual_sig_units = nan( size(unit_labels, 1), 1 );
m2_mutual_sig_units = nan( size(unit_labels, 1), 1 );

[I, C] = findall( psth_labels, {'uuid'}, mask );
fprintf('Method used for stat comparison: %s\n', method);
for i = 1:numel(I)
    inds = I{i};
    roi = 'eyes_nf';
    m1_excl_inds = find( psth_labels, {'m1', roi}, inds );
    m1_init_inds = find( psth_labels, {'mutual', 'm1_initiated', roi}, inds );
    m2_excl_inds = find( psth_labels, {'m2', roi}, inds );
    m2_init_inds = find( psth_labels, {'mutual', 'm2_initiated', roi}, inds );
    m1_mutual_sig_units(i) = eisg.mean_fr.stat_compare_mean_fr(...
        psth_matrix, t, m1_excl_inds, m1_init_inds, method, analysis_time_window );
    m2_mutual_sig_units(i) = eisg.mean_fr.stat_compare_mean_fr(...
        psth_matrix, t, m2_excl_inds, m2_init_inds, method, analysis_time_window );
end

%%

sig_units = m2_mutual_sig_units;
sig_units(isnan( sig_units )) = 0;
use_labels = unit_labels;
mask = rowmask( use_labels );

[~, cell_type_I, cell_types] = keepeach(use_labels', {'cell-type', 'region'});
cell_type_count = cellfun(@(x) numel(sig_units(x)), cell_type_I);
sig_cell_count = cellfun(@(x) sum( sig_units(x) == 1 ), cell_type_I);
sig_cell_frac = cellfun(@(x) sum( sig_units(x) == 1 )/numel(x), cell_type_I);
prop_table = table(cell_type_count, sig_cell_count, sig_cell_frac, 'RowNames', fcat.strjoin(cell_types, '_'));
disp('Proportion table for mutual gaze encoding:');
fprintf('Time window for mean: %0.2fs to %0.2fs\n', analysis_time_window(1), analysis_time_window(2) );
disp(prop_table);

[tbl, reg_labels] = do_prop_tests(sig_units, use_labels', mask);
disp('Chi-Square Test comparing celltype proportions:');
disp(tbl);

%% Helper Functions
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

