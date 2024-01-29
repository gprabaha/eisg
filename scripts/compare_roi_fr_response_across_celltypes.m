% Workspace Cleanup
clc;
clear;

%% Script Parameters
do_psth_extraction = false;
smoothen_psth = false;

% Validity filter for units
validity_filter = {'valid-unit', 'maybe-valid-unit'};

% For analysis
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

%% Declare PSTH Extraction Parameters
min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;
pre_time_range = [-0.5 0]; % in seconds
post_time_range = [0 0.5];

%% Compute/Load PSTH and Associated Parameters
if do_psth_extraction
    disp('Computing PSTH...');
    rois = {'face', 'eyes_nf', 'whole_face', 'right_nonsocial_object_whole_face_matched'};
    evt_mask = find(events.labels, [{'m1'}, rois]);
    spk_mask = find(spike_labels, validity_filter); 
    [psth_matrix, psth_labels, t] = eisg.psth.compute_psth(...
        unit_spike_ts, spike_labels, spk_mask, ...
        evts, events.labels, evt_mask, ...
        min_t, max_t, bin_width);
    pre_time = t > pre_time_range(1) & t <= pre_time_range(2);
    post_time = t > post_time_range(1) & t <= post_time_range(2);
    % Save variables
    disp('Saving PSTH...');
    save(fullfile(data_p, 'binned_unit_psth_social_gaze.mat')...
        , 'psth_matrix', 'psth_labels', 't'...
        );
else
    disp('Loading saved PSTH...');
    loaded_data = load(fullfile(data_p, 'binned_unit_psth_social_gaze.mat'));
    psth_matrix = loaded_data.psth_matrix;
    psth_labels = loaded_data.psth_labels;
    t = loaded_data.t;
    pre_time = t > pre_time_range(1) & t <= pre_time_range(2);
    post_time = t > post_time_range(1) & t <= post_time_range(2);
end

%% Smoothen Data
% !!
% !! Smoothening data or not has a impact on the stats
% !! ACC vs BLA prop test works out only if unsmoothened
% !!
if smoothen_psth
    disp('Smoothening PSTH...'); 
    psth_matrix = eisg.psth.smoothen_psth(psth_matrix);
else
    disp('Using raw PSTH without smoothening...');
end

%% Calculate Mean FR per Unit
disp('Calculating mean pre- and post-fr per trial...');
spike_pre = nanmean(psth_matrix(:, pre_time), 2);
spike_post = nanmean(psth_matrix(:, post_time), 2);

%% ANOVA Tests
disp('Performing ANOVA tests...')
mask = find(psth_labels, {'eyes_nf', 'face', 'right_nonsocial_object_whole_face_matched'});
[anova_labels, anova_I] = keepeach(psth_labels', 'uuid', mask);
anova_outs = eisg.anova.hanova(spike_pre, psth_labels', mask, anova_I);
anova_outs_post = eisg.anova.hanova(spike_post, psth_labels', mask, anova_I);
anova_ps = cate1(anova_outs.ps);
anova_post_ps = cate1(anova_outs_post.ps);
disp('Done');

%% Cell Type Analysis - ANOVA Factor 2 (Face vs Obj)
disp('Printing results...');
factor_ind = 2;
use_labels_face_obj = anova_labels';
sig_cells_face_obj = anova_ps(:, factor_ind) < 0.05 | anova_post_ps(:, factor_ind) < 0.05;


% Find indices of categories to exclude
chi2_mask = findnone(anova_labels, excluded_categories);

% Apply mask to anova_labels
[~, cell_type_I_face_obj, cell_types_face_obj] = keepeach(anova_labels'...
    , {'cell-type', 'region'}, chi2_mask);

cell_type_count_face_obj = cellfun(@(x) numel(sig_cells_face_obj(x)), cell_type_I_face_obj);
sig_cell_count_face_obj = cellfun(@(x) nnz(sig_cells_face_obj(x)), cell_type_I_face_obj);
cell_type_freqs_face_obj = cellfun(@(x) nnz(sig_cells_face_obj(x))/numel(x), cell_type_I_face_obj);

prop_table_face_obj = table(cell_type_count_face_obj...
    , sig_cell_count_face_obj, cell_type_freqs_face_obj, 'RowNames'...
    , fcat.strjoin(cell_types_face_obj, '_'));
disp('Proportion Table for Factor 2 (Face vs Obj):');
disp(prop_table_face_obj);
[tbl_face_obj, reg_labels_face_obj] = do_prop_tests(sig_cells_face_obj...
    , use_labels_face_obj', findnone(use_labels_face_obj...
    , {'outlier', 'ofc', 'dmpfc'}));
disp('Chi-Square Test for Factor 2 (Face vs Obj):');
disp(tbl_face_obj);

plot_cell_type_pie_by_region(sig_cells_face_obj, cell_types_face_obj, cell_type_I_face_obj, 'Face vs Obj');

%% Pie charts

% Add scripts here to convert the tables to piecharts

%% Cell Type Analysis - ANOVA Factor 1 (Eye vs Non-eye Face)
factor_ind = 1;
use_labels_eye_nef = fcat.from(anova_labels(find(sig_cells_face_obj), :), getcats(anova_labels))';
sig_cells_eye_nef = sig_cells_face_obj & (anova_ps(:, factor_ind) < 0.05 | anova_post_ps(:, factor_ind) < 0.05);
sig_cells_eye_nef = sig_cells_eye_nef(sig_cells_face_obj);

chi2_mask = findnone(use_labels_eye_nef, excluded_categories);
[~, cell_type_I_eye_nef, cell_types_eye_nef] = keepeach(use_labels_eye_nef'...
    , {'cell-type', 'region'}, chi2_mask);

cell_type_count_eye_nef = cellfun(@(x) numel(sig_cells_eye_nef(x)), cell_type_I_eye_nef);
sig_cell_count_eye_nef = cellfun(@(x) nnz(sig_cells_eye_nef(x)), cell_type_I_eye_nef);
cell_type_freqs_eye_nef = cellfun(@(x) nnz(sig_cells_eye_nef(x))/numel(x), cell_type_I_eye_nef);

prop_table_eye_nef = table(cell_type_count_eye_nef, sig_cell_count_eye_nef...
    , cell_type_freqs_eye_nef, 'RowNames', fcat.strjoin(cell_types_eye_nef, '_'));
disp('Proportion Table for Factor 1 (Eye vs Non-eye Face):');
disp(prop_table_eye_nef);
[tbl_eye_nef, reg_labels_eye_nef] = do_prop_tests(sig_cells_eye_nef...
    , use_labels_eye_nef', chi2_mask);
disp('Chi-Square Test for Factor 1 (Eye vs Non-eye Face):');
disp(tbl_eye_nef);

plot_cell_type_pie_by_region(sig_cells_face_obj, cell_types_face_obj, cell_type_I_face_obj, 'Eye vs Non-eye Face');

%% Ranksum Tests
disp('Performing ranksum tests...');
pre_rs_outs = dsp3.ranksum( spike_pre, psth_labels', {'uuid'}, 'whole_face', 'right_nonsocial_object_whole_face_matched' );
post_rs_outs = dsp3.ranksum( spike_post, psth_labels', {'uuid'}, 'whole_face', 'right_nonsocial_object_whole_face_matched' );

pre_rs_ps = cellfun(@(x) x.p, pre_rs_outs.rs_tables);
post_rs_ps = cellfun(@(x) x.p, post_rs_outs.rs_tables);

pre_rs_zs = cellfun(@(x) x.zval, pre_rs_outs.rs_tables);
post_rs_zs = cellfun(@(x) x.zval, post_rs_outs.rs_tables);

%% Proportion Comparison Across Cell-types (Face vs Obj) for Ranksum
disp('Printing results...');
use_labels_rs = pre_rs_outs.rs_labels';
sig_cells = pre_rs_ps < 0.05 | post_rs_ps < 0.05;
chi2_mask = findnone( use_labels_rs', excluded_categories );
[~, cell_type_I_rs, cell_types_rs] = keepeach(use_labels_rs', {'cell-type', 'region'}, chi2_mask);
cell_type_count_rs = cellfun(@(x) numel(sig_cells(x)), cell_type_I_rs);
sig_cell_count_rs = cellfun(@(x) nnz(sig_cells(x)), cell_type_I_rs);
cell_type_freqs_rs = cellfun(@(x) nnz(sig_cells(x))/numel(x), cell_type_I_rs);
prop_table_rs = table(cell_type_count_rs, sig_cell_count_rs, cell_type_freqs_rs, 'RowNames', fcat.strjoin(cell_types_rs, '_'));
disp('Proportion Table for Ranksum Tests:');
disp(prop_table_rs);

[tbl_rs, reg_labels_rs] = do_prop_tests(sig_cells, use_labels_rs', chi2_mask);
disp('Chi-Square Test for Face vs Obj Ranksum Tests:');
disp(tbl_rs);

%%
function plot_cell_type_pie_by_region(sig_cells, cell_types_face_obj, cell_type_I_face_obj, additional_title)
    % Get unique regions
    regions = unique(cell_types_face_obj(2, :));

    % Loop through each region
    for r = 1:numel(regions)
        current_region = regions{r};
        
        % Get indices for current region
        region_indices = find(strcmp(cell_types_face_obj(2, :), current_region));

        % Create a new figure for the current region
        figure;
        sgtitle(['Significant Cell Fraction by Cell Type in ' current_region ' - ' additional_title]);

        % Loop through each cell type within the region
        for i = 1:numel(region_indices)
            % Get cell type
            cell_type = cell_types_face_obj{1, region_indices(i)};
            
            % Get indices for current cell type
            current_indices = cell_type_I_face_obj{region_indices(i)};

            % Calculate fraction of significant cells for current cell type within the region
            num_sig_cells = sum(sig_cells(current_indices));
            num_total_cells = numel(current_indices);
            fraction_sig_cells = num_sig_cells / num_total_cells;

            % Format percentage
            percentage_sig = sprintf('Percentage significant: %.2f%% (%d/%d units)', fraction_sig_cells * 100, num_sig_cells, num_total_cells);

            % Plot pie chart in a subplot
            subplot(1, 2, i);
            pie([fraction_sig_cells, 1 - fraction_sig_cells], {'Significant', 'Not Significant'});
            title({cell_type; percentage_sig});
        end
    end
end




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
