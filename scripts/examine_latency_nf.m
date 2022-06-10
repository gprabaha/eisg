data_p = '~/Downloads';

%%
sorted = shared_utils.io.fload( fullfile(data_p, 'sorted_neural_data_social_gaze.mat') );
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );
ct_labels = shared_utils.io.fload( fullfile(data_p, 'celltype_labels_combined.mat') );
%%
[unit_spike_ts, unit_wfs, spike_labels] = linearize_sorted( sorted );
spike_labels = apply_cell_type_labels(...
  spike_labels, ct_labels.label_mat, ct_labels.label_mat_cols );

events = add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

%%  psth

rois = {'whole_face', 'right_nonsocial_object_whole_face_matched'};

evt_mask = find( events.labels, [{'m1'}, rois] );
spk_mask = find( spike_labels, {'valid-unit', 'maybe-valid-unit'} );

min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;

[psth_matrix, psth_labels, t] = compute_psth(...
  unit_spike_ts, spike_labels, spk_mask, evts, events.labels, evt_mask, min_t, max_t, bin_width );

method = 'square_wave';
psth_matrix_smoothened = eisg.util.get_overlapping_50ms_psth_from_nonoverlapping_10ms_psth( psth_matrix, method );

%%  baseline psth stats

% Because the baseline psth matrix would be too big to fit in memory,
% process cells from one session at a time, keeping only the mean and sem
% of spike counts from each cell.

% base_rois = { 'everywhere', 'right_nonsocial_object_whole_face_matched', 'whole_face' };
% base_rois = { 'right_nonsocial_object_whole_face_matched', 'whole_face' };
base_rois = { 'everywhere' };

evt_mask = find( events.labels, [{'m1'}, base_rois] );
spk_mask = find( spike_labels, {'valid-unit', 'maybe-valid-unit'} );

if ( 0 )  % debug - process only one session
  evt_mask = find( events.labels, '01082019', evt_mask );
end
% 
% [base_stats, base_stat_labs] = baseline_psth_stats( ...
%   evts, events.labels, evt_mask, spike_ts, spike_labels, spk_mask ...
%   , 'min_t', min_t ...
%   , 'max_t', max_t ...
%   , 'bin_width', bin_width ...
% );

[base_stats2, base_stat_labs] = baseline_psth_stats( ...
  evts, events.labels, evt_mask, unit_spike_ts, spike_labels, spk_mask ...
  , 'min_t', min_t ...
  , 'max_t', max_t ...
  , 'bin_width', bin_width ...
);

%%

pl = plotlabeled.make_common();
base_stat_labs = apply_cell_type_labels(base_stat_labs, ct_labels.label_mat, ct_labels.label_mat_cols);
% axs = pl.violinplot( base_stats(:, 1), base_stat_labs, {'cell-type'}, 'region' );
axs = pl.violinplot( base_stats2(:, 1), base_stat_labs, {'cell-type'}, 'region' );

%%  mean psth

[mean_labs, mean_I] = retaineach( psth_labels, {'uuid', 'looks_by', 'roi'} );
psth_means = bfw.row_nanmean( psth_matrix_smoothened, mean_I );

%%  peak times from rows

n_devs = 4; % detect peak at first time above mean + std * n_devs
peak_indices = find_first( above_sem(psth_means, n_devs) );
peaks = to_peak_times( peak_indices, t );

%%  peak times from baseline

n_devs = 4; % detect peak at first time above mean + std * n_devs
peak_I = bfw.find_combinations( base_stat_labs, mean_labs(:, 'uuid')' );
above_base = above_sems( psth_means, peak_I, base_stats(:, 1), base_stats(:, 2), n_devs );
peak_indices = find_first( above_base );
peaks = to_peak_times( peak_indices, t );

%%  plot cumuluative histogram of latencies

save_p = fullfile( eisg.util.project_path, 'data/plots/latency/cum_hist', dsp3.datedir );
do_save = false;

assert_ispair( peak_indices, mean_labs );
peak_mat = to_peak_matrix( peak_indices, numel(t) );
[lat_labs, lat_I] = retaineach( mean_labs, {'region', 'looks_by', 'cell-type', 'roi'} );
peak_hist = to_histogram( peak_mat, lat_I );

plt_mask = pipe( rowmask(lat_labs) ...
  , @(m) findnone(lat_labs, '<cell-type>', m) ...
);

fig_I = findall( lat_labs, 'region', plt_mask );

for i = 1:numel(fig_I)
  fi = fig_I{i};
  
  pl = plotlabeled.make_common();
  pl.x = t;
  pcats = {'cell-type', 'region', 'looks_by'};
  gcats = {'roi'};
  axs = pl.lines( peak_hist(fi, :), lat_labs(fi), gcats, pcats );
  shared_utils.plot.set_ylims( axs, [0, 1] );
  
  if ( do_save )
    dsp3.req_savefig( gcf, save_p, prune(lat_labs(fi)), pcats );
  end
end

%%  plot histogram of latencies

save_p = fullfile( eisg.util.project_path, 'data/plots/latency/hist', dsp3.datedir );
do_save = true;

plt_mask = pipe( rowmask(mean_labs) ...
  , @(m) intersect(m, find(~isnan(peaks))) ...
  , @(m) findnone(mean_labs, '<cell-type>', m) ...
);

fig_I = findall( mean_labs, {'region'}, plt_mask );
ct_labs = combs( mean_labs, 'cell-type', plt_mask );

for i = 1:numel(fig_I)
  fi = fig_I{i};
  
  pl = plotlabeled.make_common();
  pl.hist_add_summary_line = true;
  pl.summary_func = @nanmedian;
  pl.panel_order = ct_labs;
  pl.shape = [3, 2];
  
  pcats = {'roi', 'cell-type', 'region'};  
  plt_labs = mean_labs(fi);
  plt_dat = peaks(fi);
  axs = pl.hist( plt_dat, plt_labs, pcats, 20 );
  
  if ( do_save )
    dsp3.req_savefig( gcf, save_p, prune(plt_labs), pcats );
  end
end

%%

function peaks = to_peak_times(peak_indices, t)

peaks = peak_indices;
peaks(peaks == 0) = nan;
peaks(~isnan(peaks)) = t(peaks(~isnan(peaks)));

end

function [base_stats, base_stat_labs] = baseline_psth_stats(...
  evts, evt_labels, evt_mask, spike_ts, spike_labels, spk_mask, varargin)

defaults = struct();
defaults.min_t = -0.5;
defaults.max_t = 0.5;
defaults.bin_width = 0.05;

params = shared_utils.general.parsestruct( defaults, varargin );

base_stats = [];
base_stat_labs = fcat();

[evt_I, evt_C] = findall( evt_labels, 'session', evt_mask );
for i = 1:numel(evt_I)
  fprintf( '\n Session %s (%d of %d)', evt_C{i}, i, numel(evt_I) );
  
  sesh_evt_mask = evt_I{i};
  sesh_spk_mask = find( spike_labels, evt_C(:, i), spk_mask );
  spk_I = findall( spike_labels, 'uuid', sesh_spk_mask );
  spk_Is = shared_utils.vector.bin( spk_I, 4 );
  
  for j = 1:numel(spk_Is)
    fprintf( '\n\t Group (%d of %d)', j, numel(spk_Is) );
    
    si = vertcat( spk_Is{j}{:} );
    % changes made here
    [base_psth_matrix, base_psth_labels, t] = compute_psth(...
      spike_ts, spike_labels, si, evts, evt_labels, sesh_evt_mask ...
      , params.min_t, params.max_t, params.bin_width );
    time_snippet = t >= 0.05 & t < 0.4;
    base_psth_matrix = mean( base_psth_matrix( :, time_snippet), 2 );
    
    [unit_labs, unit_I] = retaineach( base_psth_labels, {'uuid', 'looks_by'} );

  %   sigmas = cellfun( @(x) std(columnize(base_psth_matrix(x))), unit_I );
    sigmas = cellfun( @(x) sem_all(columnize(base_psth_matrix(x))), unit_I );
    means = cellfun( @(x) mean(columnize(base_psth_matrix(x))), unit_I );
    base_stats = [ base_stats; [means, sigmas] ];

    append( base_stat_labs, unit_labs );
  end
end

end