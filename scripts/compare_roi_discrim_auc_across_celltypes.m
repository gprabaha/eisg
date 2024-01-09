%% Workspace Cleanup
clc;
clear;

%% Script Parameters
do_psth_extraction = false;
smoothen_psth = true;

% Validity filter for units
validity_filter = {'valid-unit', 'maybe-valid-unit'};

% For plotting
excluded_categories = {'outlier', 'ofc', 'dmpfc'};

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

%% Face vs Obj AUC Calculation
roi_a = 'whole_face';
roi_b = 'right_nonsocial_object_whole_face_matched';
[auc_f_o, z_scored_auc_f_o, auc_labels_f_o] = eisg.auc.calculate_roi_comparison_auc(...
    psth_matrix, psth_labels, roi_a, roi_b...
    );

%% Plot Face vs Obj AUC Timeseries for ACCg
region = 'acc';
eisg.plot.plot_zscored_auc_timeseries_across_celltypes(...
    t, z_scored_auc_f_o, auc_labels_f_o, region, excluded_categories);
sgtitle(['Face vs Obj AUC Timeseries for: ' region]);

%% Plot Face vs Obj AUC Timeseries for BLA
region = 'bla';
eisg.plot.plot_zscored_auc_timeseries_across_celltypes(...
    t, z_scored_auc_f_o, auc_labels_f_o, region, excluded_categories);
sgtitle(['Face vs Obj AUC Timeseries for: ' region]);

%% Eyes vs Non-eye Face AUC Calculation
roi_a = 'eyes_nf';
roi_b = 'face';
[auc_enf_f, z_scored_auc_enf_f, auc_labels_enf_f] = eisg.auc.calculate_roi_comparison_auc(...
    psth_matrix, psth_labels, roi_a, roi_b...
    );

%% Plot Eyes vs Non-eye Face AUC Timeseries for ACCg
region = 'acc';
eisg.plot.plot_zscored_auc_timeseries_across_celltypes(...
    t, z_scored_auc_enf_f, auc_labels_enf_f, region, excluded_categories);
sgtitle(['Eyes vs Non-eye Face AUC Timeseries for: ' region]);

%% Plot Eyes vs Non-eye Face AUC Timeseries for BLA
region = 'bla';
eisg.plot.plot_zscored_auc_timeseries_across_celltypes(...
    t, z_scored_auc_enf_f, auc_labels_enf_f, region, excluded_categories);
sgtitle(['Eyes vs Non-eye Face AUC Timeseries for: ' region]);
 

%{
%% Mean AUC/unit comparison: ACC and BLA

% Calculate mean and SEM for ACC
mean_acc_broad = nanmedian( nanmean( abs_zscore_auc_broad{2}, 2 ) );
sem_acc_broad = sem( nanmean( abs_zscore_auc_broad{2}, 2 ) );
mean_acc_narrow = nanmedian( nanmean(abs_zscore_auc_narrow{2}, 2 ) );
sem_acc_narrow = sem( nanmean( abs_zscore_auc_narrow{2}, 2 ) );

% Calculate mean and SEM for BLA
mean_bla_broad = nanmedian( nanmean(abs_zscore_auc_broad{1}, 2 ) );
sem_bla_broad = sem(nanmean( abs_zscore_auc_broad{1}, 2) );
mean_bla_narrow = nanmedian( nanmean(abs_zscore_auc_narrow{1}, 2 ) );
sem_bla_narrow = sem(nanmean( abs_zscore_auc_narrow{1}, 2) );

% Bar graph parameters
means_acc = [mean_acc_broad, mean_acc_narrow];
sems_acc = [sem_acc_broad, sem_acc_narrow];
means_bla = [mean_bla_broad, mean_bla_narrow];
sems_bla = [sem_bla_broad, sem_bla_narrow];
xticklabels = {'ACC', 'BLA'};

% Create figure
figure;

% Plot means for ACC
bar(1, means_acc(1), 'BarWidth', 0.5, 'FaceColor', 'r');
hold on;
bar(2, means_acc(2), 'BarWidth', 0.5, 'FaceColor', 'b');

% Plot means for BLA
bar(4, means_bla(1), 'BarWidth', 0.5, 'FaceColor', 'r');
bar(5, means_bla(2), 'BarWidth', 0.5, 'FaceColor', 'b');
hold off;

% Plot standard errors for ACC
hold on;
errorbar([1, 2], means_acc, sems_acc, 'k.', 'LineWidth', 1, 'CapSize', 10);

% Plot standard errors for BLA
errorbar([4, 5], means_bla, sems_bla, 'k.', 'LineWidth', 1, 'CapSize', 10);
hold off;

% Customize the plot
title('Mean AUC/unit Comparison - Broad vs Narrow');
xlabel('Category');
ylabel('Mean AUC');
legend('Broad', 'Narrow', 'Location', 'best');
set(gca, 'XTick', [1.5, 4.5], 'XTickLabel', xticklabels);
grid on;


% Check the stats
p_acc = ranksum( nanmedian( abs_zscore_auc_broad{2} ), nanmedian( abs_zscore_auc_narrow{2} ) );
median_acc_broad = nanmedian( nanmean( abs_zscore_auc_broad{2}, 2 ) );
median_acc_narrow = nanmedian( nanmean( abs_zscore_auc_narrow{2}, 2 ) );
fprintf('Median AUC/bin ACC broad: %0.3f; ACC narrow: %0.3f; p: %f\n', median_acc_broad, median_acc_narrow, p_acc );

p_bla = ranksum( nanmedian( abs_zscore_auc_broad{1} ), nanmedian( abs_zscore_auc_narrow{1} ) );
median_bla_broad = nanmedian( nanmean( abs_zscore_auc_broad{1}, 2 ) );
median_bla_narrow = nanmedian( nanmean( abs_zscore_auc_narrow{1}, 2 ) );
fprintf('Median AUC/bin BLA broad: %0.3f; ACC narrow: %0.3f; p: %f\n', median_bla_broad, median_bla_narrow, p_bla );


p_acc = ranksum( nanmean( abs_zscore_auc_broad{2}, 2 ), nanmean( abs_zscore_auc_narrow{2}, 2 ) );
mean_acc_broad = nanmean( nanmean( abs_zscore_auc_broad{2} ) );
mean_acc_narrow = nanmean( nanmean( abs_zscore_auc_narrow{2} ) );
fprintf('Median AUC/unit ACC broad: %0.3f; ACC narrow: %0.3f; p: %f\n', mean_acc_broad, mean_acc_narrow, p_acc );

p_bla = ranksum( nanmean( abs_zscore_auc_broad{1}, 2 ), nanmean( abs_zscore_auc_narrow{1}, 2 ) );
mean_bla_broad = nanmedian( nanmean( abs_zscore_auc_broad{1}, 2 ) );
mean_bla_narrow = nanmedian( nanmean( abs_zscore_auc_narrow{1}, 2 ) );
fprintf('Median AUC/unit BLA broad: %0.3f; ACC narrow: %0.3f; p: %f\n', mean_bla_broad, mean_bla_narrow, p_bla );


%% Max AUC/unit comparison: ACC and BLA

% Calculate max and SEM for ACC
max_acc_broad = nanmean(max(abs_zscore_auc_broad{2}, [], 2));
sem_acc_broad = sem(max(abs_zscore_auc_broad{2}, [], 2));
max_acc_narrow = nanmean(max(abs_zscore_auc_narrow{2}, [], 2));
sem_acc_narrow = sem(max(abs_zscore_auc_narrow{2}, [], 2));

% Calculate max and SEM for BLA
max_bla_broad = nanmean(max(abs_zscore_auc_broad{1}, [], 2));
sem_bla_broad = sem(max(abs_zscore_auc_broad{1}, [], 2));
max_bla_narrow = nanmean(max(abs_zscore_auc_narrow{1}, [], 2));
sem_bla_narrow = sem(max(abs_zscore_auc_narrow{1}, [], 2));

% Bar graph parameters
means_acc = [max_acc_broad, max_acc_narrow];
sems_acc = [sem_acc_broad, sem_acc_narrow];
means_bla = [max_bla_broad, max_bla_narrow];
sems_bla = [sem_bla_broad, sem_bla_narrow];
xticklabels = {'ACC', 'BLA'};

% Create figure
figure;

% Plot means for ACC
bar(1, means_acc(1), 'BarWidth', 0.5, 'FaceColor', 'r');
hold on;
bar(2, means_acc(2), 'BarWidth', 0.5, 'FaceColor', 'b');

% Plot means for BLA
bar(4, means_bla(1), 'BarWidth', 0.5, 'FaceColor', 'r');
bar(5, means_bla(2), 'BarWidth', 0.5, 'FaceColor', 'b');
hold off;

% Plot standard errors for ACC
hold on;
errorbar([1, 2], means_acc, sems_acc, 'k.', 'LineWidth', 1, 'CapSize', 10);

% Plot standard errors for BLA
errorbar([4, 5], means_bla, sems_bla, 'k.', 'LineWidth', 1, 'CapSize', 10);
hold off;

% Customize the plot
title('Max AUC Comparison - Broad vs Narrow');
xlabel('Category');
ylabel('Max AUC');
legend('Broad', 'Narrow', 'Location', 'best');
set(gca, 'XTick', [1.5, 4.5], 'XTickLabel', xticklabels);
grid on;

p_acc = ranksum( max( abs_zscore_auc_broad{2}, [], 2 ), max( abs_zscore_auc_narrow{2}, [], 2 ) );
median_max_acc_broad = nanmedian( max( abs_zscore_auc_broad{2}, [], 2 ) );
median_max_acc_narrow = nanmedian( max( abs_zscore_auc_narrow{2}, [], 2 ) );
fprintf('Median max AUC val ACC; broad: %0.3f; ACC narrow: %0.3f; p: %f\n', median_max_acc_broad, median_max_acc_narrow, p_acc );

p_bla = ranksum( max( abs_zscore_auc_broad{1}, [], 2 ), max( abs_zscore_auc_narrow{1}, [], 2 ) );
median_max_bla_broad = nanmedian( max( abs_zscore_auc_broad{1}, [], 2 ) );
median_max_bla_narrow = nanmedian( max( abs_zscore_auc_narrow{2}, [], 2 ) );
fprintf('Median max AUC val BLA; broad: %0.3f; ACC narrow: %0.3f; p: %f\n', median_max_bla_broad, median_max_bla_narrow, p_bla );

%%
function [semval] = sem(vector_data)

% Recall that s.e.m. = std(x)/sqrt(length(x));
nonan_std = nanstd(vector_data);
nonan_len = length(vector_data(~isnan(vector_data)));

% Plug in values
semval = nonan_std / sqrt(nonan_len);

end

%}