data_p = '/Users/prabaha/repositories/eisg/processed_data';

sorted = shared_utils.io.fload( fullfile(data_p, 'sorted_neural_data_social_gaze.mat') );
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );
ct_labels = shared_utils.io.fload( fullfile(data_p, 'celltype_labels_combined.mat') );

[spike_ts, spike_labels] = linearize_sorted( sorted );
spike_labels = apply_cell_type_labels(...
  spike_labels, ct_labels.label_mat, ct_labels.label_mat_cols );

events = add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

%%  psth

rois = { 'whole_face', 'right_nonsocial_object_whole_face_matched' };

evt_mask = find( events.labels, [{'m1'}, rois] );
spk_mask = find( spike_labels, {'valid-unit', 'maybe-valid-unit'} );

min_t = -0.5;
max_t = 0.5;
bin_width = 0.05;

[psth_matrix, psth_labels, t] = compute_psth(...
    spike_ts, spike_labels, spk_mask ...
  , evts, events.labels, evt_mask ...
  , min_t, max_t, bin_width );

%% compare using ANOVA instead of ranksum

% get the roi-pair indices from psth labels
rois = combs( psth_labels, 'roi' );
roi_pairs = bfw.pair_combination_indices( numel(rois) );

% get the a-s and b-s to compare
as = arrayfun( @(x) rois{x}, roi_pairs(:, 1), 'un', 0 );
bs = arrayfun( @(x) rois{x}, roi_pairs(:, 2), 'un', 0 );

% analysis for POST time window

% post time window
t_win_post = [0, 0.5];
t_mask_post = t >= t_win_post(1) & t < t_win_post(2);
t_mean_post = mean( psth_matrix(:, t_mask_post), 2 );

% ranksum tests comparing a vs b, for each a in as and b in bs, for post
% fixation epoch
[rs_tables_post, rs_labels_post] = ranksums( ...
  t_mean_post, psth_labels, {'uuid', 'looks_by'}, as, bs, 'roi' );

% extract p value from ranksum test for post fixation epoch
rs_post_ps = cellfun( @(x) x.p, rs_tables_post );
is_sig_post = rs_post_ps < 0.05;
sig_cells_post = cellstr( rs_labels_post, {'region', 'uuid'}, find(is_sig_post) );

% analysis for PRE time window

% pre time window
t_win_pre = [-0.5, 0];
t_mask_pre = t >= t_win_pre(1) & t < t_win_pre(2);
t_mean_pre = mean( psth_matrix(:, t_mask_pre), 2 );

% ranksum tests comparing a vs b, for each a in as and b in bs, for pre
% fixation epoch
[rs_tables_pre, rs_labels_pre] = ranksums( ...
  t_mean_pre, psth_labels, {'uuid', 'looks_by'}, as, bs, 'roi' );

% extract p value from ranksum test for post fixation epoch
rs_pre_ps = cellfun( @(x) x.p, rs_tables_pre );
is_sig_pre = rs_pre_ps < 0.05;
sig_cells_pre = cellstr( rs_labels_pre, {'region', 'uuid'}, find(is_sig_pre) );

% PRE and POST combined
is_sig = is_sig_pre | is_sig_post;
sig_cells = cellstr( rs_labels_post, {'region', 'uuid'}, find(is_sig) );

% proportions of significant cells per region and roi
[prop_labs, prop_I] = retaineach( rs_labels_post, {'region', 'cell-type', 'looks_by', 'roi'} );
props = proportions( is_sig, prop_I );

plt_mask = pipe( rowmask(prop_labs) ...
  , @(m) findnone(prop_labs, '<cell-type>', m) ... % remove unclassified celltypes
);

pl = plotlabeled.make_common();
axs = pl.bar( props(plt_mask), prop_labs(plt_mask), {}, 'cell-type', {'roi', 'looks_by', 'region'} );
ylabel( axs(1), '% Significant a vs b.' );

%% update face-obj discriminating neurons in spike_labels

spike_labels = apply_face_obj_sig_labels(...
  spike_labels, sig_cells );

% %%  mean psth + peak times
% 
% % This has to be SEMs, not STDs
% n_devs = 1.5; % detect peak at first time above mean + std * n_devs
% 
% [mean_labs, mean_I] = retaineach( psth_labels, {'uuid', 'looks_by', 'roi'} );
% psth_means = bfw.row_nanmean( psth_matrix, mean_I );
% 
% peaks = find_first( above_std(psth_means, n_devs) );
% peaks(peaks == 0) = nan;
% peaks(~isnan(peaks)) = t(peaks(~isnan(peaks)));
% 
% 
% %% plot histogram of latencies
% 
% save_p = fullfile( eisg.util.project_path, 'data/plots/latency/hist', dsp3.datedir );
% do_save = true;
% 
% plt_mask = pipe( rowmask(mean_labs) ...
%   , @(m) intersect(m, find(~isnan(peaks))) ...
%   , @(m) findnone(mean_labs, '<cell-type>', m) ... % remove unclassified celltypes
%   , @(m) findnone(mean_labs, '<face-vs-obj>', m) ... % remove units that are not face-obj discriminating
% );
% 
% fig_I = findall( mean_labs, {'region'}, plt_mask );
% ct_labs = combs( mean_labs, 'cell-type', plt_mask );
% 
% for i = 1:numel(fig_I)
%   fi = fig_I{i};
%   
%   pl = plotlabeled.make_common();
%   pl.hist_add_summary_line = true;
%   pl.summary_func = @nanmedian;
%   pl.panel_order = ct_labs;
%   pl.shape = [3, 2];
%   
%   pcats = {'roi', 'cell-type', 'region'};  
%   plt_labs = mean_labs(fi);
%   plt_dat = peaks(fi);
%   axs = pl.hist( plt_dat, plt_labs, pcats, 20 );
%   
%   if ( do_save )
%     dsp3.req_savefig( gcf, save_p, prune(plt_labs), pcats );
%   end
% end
% 
% %%  discriminate between pairs of rois
% 
% t_win_post = [0, 0.5];
% t_mask_post = t >= t_win_post(1) & t < t_win_post(2);
% t_mean_post = mean( psth_matrix(:, t_mask_post), 2 );
% 
% rois = combs( psth_labels, 'roi' );
% roi_pairs = bfw.pair_combination_indices( numel(rois) );
% 
% as = arrayfun( @(x) rois{x}, roi_pairs(:, 1), 'un', 0 );
% bs = arrayfun( @(x) rois{x}, roi_pairs(:, 2), 'un', 0 );
% 
% % discriminate between a and b, for each a in as and b in bs
% each_I = findall( psth_labels, {'uuid', 'looks_by'} );
% [fit_res, fit_labels] = binary_fitcdiscrs( ...
%   t_mean_post, psth_labels, each_I, as, bs, 'roi' );
% 
% %%  plot results
% 
% % extract p value from classification
% ps = cellfun( @(x) x.p, fit_res );
% is_sig_post = ps < 0.05;
% 
% % proportions of significant cells per region and roi
% mask = findnone( fit_labels, {'<cell-type>'} );
% [prop_labs, prop_I] = retaineach(...
%   fit_labels, {'region', 'looks_by', 'roi', 'cell-type'}, mask );
% props = proportions( is_sig_post, prop_I );
% 
% pl = plotlabeled.make_common();
% axs = pl.bar( props, prop_labs, {}, 'region', {'roi', 'looks_by', 'cell-type'} );
% ylabel( axs(1), 'Prop Significant a vs b.' );
