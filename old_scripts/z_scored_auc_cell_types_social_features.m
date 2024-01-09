

clc;
clear;
data_p = fullfile( eisg.util.project_path, 'processed_data');
sorted = shared_utils.io.fload( fullfile(data_p,...
  'sorted_neural_data_social_gaze.mat') );
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );
ct_labels = load_cell_type_labels( data_p );

[unit_spike_ts, unit_wfs, spike_labels] = linearize_sorted( sorted );

bfw.add_monk_labels( spike_labels );

[uuid_I, uuids] = findall( spike_labels, 'uuid',...
  find(spike_labels, {'valid-unit', 'maybe-valid-unit'}) );
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

events = add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

%%

rois = { 'whole_face', 'right_nonsocial_object_whole_face_matched',...
  'eyes_nf', 'face' };
evt_mask = find( events.labels, [{'m1'}, rois] );
spk_mask = find( spike_labels, {'valid-unit', 'maybe-valid-unit'} );
min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;
[psth_matrix, psth_labels, t] = compute_psth(...
    unit_spike_ts, spike_labels, spk_mask ...
  , evts, events.labels, evt_mask ...
  , min_t, max_t, bin_width );

psth_matrix = smoothdata( psth_matrix, 2, 'movmean', 10 );

%%

[auc_labels, unit_I] = keepeach( psth_labels', 'uuid' );
roi_a = 'eyes_nf';
roi_b = 'face';
aucs = nan( numel(unit_I), size(psth_matrix, 2) );
null_aucs = cell( numel(unit_I), 1 );
perm_p = nan( size(aucs) );
z_scored_aucs = nan( size(aucs) );
parfor i = 1:numel(unit_I)
  fprintf( '%d of %d\n', i, numel(unit_I) );
  ind_a = find( psth_labels, roi_a, unit_I{i} );
  ind_b = find( psth_labels, roi_b, unit_I{i} );
  aucs(i, :) = auc_over_time( psth_matrix, ind_a, ind_b ); 
  null_aucs{i} = auc_perm_test( psth_matrix, ind_a, ind_b, 100 );
  z_scored_aucs(i, :) = (aucs(i, :) - mean( null_aucs{i}, 1 )) ./ std( null_aucs{i}, [], 1 );  
end


%% Broad

unit_mask = findnone( auc_labels, {'outlier', 'narrow', 'ofc', 'dmpfc'} );
[I, C] = findall( auc_labels, {'cell-type', 'region'}, unit_mask );
figure();
axs = plots.panels( numel(I) );
for i = 1:numel(axs)
  ind = I{i};
  abs_zscore_auc_broad{i} = abs(z_scored_aucs(ind,:));
  [max_aucs{i}, max_auc_loc{i}] = max( abs_zscore_auc_broad{i}, [], 2);
  [~,sorted_ind{i}] = sort(max_auc_loc{i});
  imagesc( axs(i), t, 1:numel(ind), abs_zscore_auc_broad{i}(sorted_ind{i}, :) );
  colorbar( axs(i) );
  set( axs(i), 'clim', [0, max(max(abs(z_scored_aucs)))] );
  red_map = [linspace(1, 1, 100)', linspace(1, 0, 100)', linspace(1, 0, 100)']; % Red
  custom_broad_map = [linspace(1, 0.7804, 100)', linspace(1, 0.3216, 100)', linspace(1, 0.1647, 100)']; % Kinda orange #c7522a
  colormap( axs(i), custom_broad_map );
  title( axs(i), strrep(fcat.strjoin(C(:, i), ' | '), '_', ' ') );
end
set(gcf, 'Position',  [200, 200, 900, 500]);

%% Narrow

unit_mask = findnone( auc_labels, {'outlier', 'broad', 'ofc', 'dmpfc'} );
[I, C] = findall( auc_labels, {'cell-type', 'region'}, unit_mask );
figure();
axs = plots.panels( numel(I) );
for i = 1:numel(axs)
  ind = I{i};
  abs_zscore_auc_narrow{i} = abs(z_scored_aucs(ind,:));
  [max_aucs{i}, max_auc_loc{i}] = max( abs_zscore_auc_narrow{i}, [], 2);
  [~,sorted_ind{i}] = sort(max_auc_loc{i});
  imagesc( axs(i), t, 1:numel(ind), abs_zscore_auc_narrow{i}(sorted_ind{i}, :) );
  colorbar( axs(i) );
  set( axs(i), 'clim', [0, max(max(abs(z_scored_aucs)))] );
  blue_map = [linspace(1, 0, 100)', linspace(1, 0, 100)', linspace(1, 1, 100)']; % Blue
  custom_narrow_map = [linspace(1, 0, 100)', linspace(1, 0.3294, 100)', linspace(1, 0.5216, 100)']; % Kinda blue #005485
  colormap( axs(i), custom_narrow_map );
  title( axs(i), strrep(fcat.strjoin(C(:, i), ' | '), '_', ' ') );
end
set(gcf, 'Position',  [200, 200, 900, 500]);


%% Mean AUC/unit comparison: ACC and BLA

% Calculate mean and SEM for ACC
mean_acc_broad = nanmean(nanmean(abs_zscore_auc_broad{2}));
sem_acc_broad = sem(nanmean(abs_zscore_auc_broad{2}, 2));
mean_acc_narrow = nanmean(nanmean(abs_zscore_auc_narrow{2}));
sem_acc_narrow = sem(nanmean(abs_zscore_auc_narrow{2}, 2));

% Calculate mean and SEM for BLA
mean_bla_broad = nanmean(nanmean(abs_zscore_auc_broad{1}));
sem_bla_broad = sem(nanmean(abs_zscore_auc_broad{1}, 2));
mean_bla_narrow = nanmean(nanmean(abs_zscore_auc_narrow{1}));
sem_bla_narrow = sem(nanmean(abs_zscore_auc_narrow{1}, 2));

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
median_acc_broad = nanmedian( nanmedian( abs_zscore_auc_broad{2} ) );
median_acc_narrow = nanmedian( nanmedian( abs_zscore_auc_narrow{2} ) );
fprintf('Median AUC/bin ACC broad: %0.3f; ACC narrow: %0.3f; p: %f\n', median_acc_broad, median_acc_narrow, p_acc );

p_bla = ranksum( nanmedian( abs_zscore_auc_broad{1} ), nanmedian( abs_zscore_auc_narrow{1} ) );
median_bla_broad = nanmedian( nanmedian( abs_zscore_auc_broad{1} ) );
median_bla_narrow = nanmedian( nanmedian( abs_zscore_auc_narrow{1} ) );
fprintf('Median AUC/bin BLA broad: %0.3f; ACC narrow: %0.3f; p: %f\n', median_bla_broad, median_bla_narrow, p_bla );


[~, p_acc] = ttest2( nanmean( abs_zscore_auc_broad{2}, 2 ), nanmean( abs_zscore_auc_narrow{2}, 2 ) );
mean_acc_broad = nanmean( nanmean( abs_zscore_auc_broad{2} ) );
mean_acc_narrow = nanmean( nanmean( abs_zscore_auc_narrow{2} ) );
fprintf('Mean AUC/unit ACC broad: %0.3f; ACC narrow: %0.3f; p: %f\n', mean_acc_broad, mean_acc_narrow, p_acc );

[~, p_bla] = ttest2( nanmean( abs_zscore_auc_broad{1}, 2 ), nanmean( abs_zscore_auc_narrow{1}, 2 ) );
mean_bla_broad = nanmean( nanmean( abs_zscore_auc_broad{1} ) );
mean_bla_narrow = nanmean( nanmean( abs_zscore_auc_narrow{1} ) );
fprintf('Mean AUC/unit BLA broad: %0.3f; ACC narrow: %0.3f; p: %f\n', mean_bla_broad, mean_bla_narrow, p_bla );


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


%% FUNCTIONS

% AUC Over Time
function aucs = auc_over_time(spikes, ind_a, ind_b)

aucs = nan( 1, size(spikes, 2) );
for i = 1:size(spikes, 2)
  spks_a = spikes(ind_a, i);
  spks_b = spikes(ind_b, i);
  auc = score_auc( spks_a, spks_b );
  aucs(i) = auc;
end

end

% AUC Perm Test
function null_aucs = auc_perm_test(spikes, ind_a, ind_b, perm_iters)

null_aucs = cell( perm_iters, 1 );
for i = 1:perm_iters
  [ia, ib] = shuffle2( ind_a, ind_b );
  null_aucs{i} = auc_over_time( spikes, ia, ib );
end
null_aucs = vertcat( null_aucs{:} );
null_aucs = sort( null_aucs );

end

% 2D Shuffle
function [ic, id] = shuffle2(ia, ib)

i = [ia; ib];
i = i(randperm(numel(i)));
ic = i(1:numel(ia));
id = i(numel(ia)+1:end);
assert( numel(ic) == numel(ia) && numel(id) == numel(ib) );

end

% AUC Score
function auc = score_auc(a, b)

t = false( numel(a) + numel(b), 1 );
t(1:numel(a)) = true;
y = [ a; b ];
auc = scoreAUC( t, y );

end

%%Fast AUC
%%Copied over from: https://www.mathworks.com/matlabcentral/fileexchange/50962-fast-auc

function auc = scoreAUC(labels,scores)
% Calculates the AUC - area under the curve.
%
% Besides being the area under the ROC curve, AUC is has a slightly 
% less known interpretation:
% If you choose a random pair of samples which is one positive and one
% negative - AUC is the probabilty that the positive-sample score is above
% the negative-sample score.
% Here we compute the AUC by counting these pairs.
% 
% auc = scoreAUC(labels,scores)
% N x 1 boolean labels
% N x 1 scores
%
% http://www.springerlink.com/content/nn141j42838n7u21/fulltext.pdf
%
% ==== Lior Kirsch 2014 ====

assert( islogical(labels) ,'labels input should be logical');
assert( isequal(size(labels),size(scores)) ,'labels and scores should have the same size');

[n,m] = size(labels);
assert( m==1, 'should have sizse of (n,1)');

num_pos = sum(labels);
num_neg = n - num_pos;
assert( ~(num_pos==0), 'no positive labels entered');
assert( ~(num_pos==n), 'no negative labels entered');


ranks = tiedrank(scores);
auc = ( sum( ranks(labels) ) - num_pos * (num_pos+1)/2) / ( num_pos * num_neg);

end

%%
function [semval] = sem(vector_data)

% Recall that s.e.m. = std(x)/sqrt(length(x));
nonan_std = nanstd(vector_data);
nonan_len = length(vector_data(~isnan(vector_data)));

% Plug in values
semval = nonan_std / sqrt(nonan_len);

end