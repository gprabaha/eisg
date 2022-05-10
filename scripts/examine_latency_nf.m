data_p = '~/Downloads';

sorted = shared_utils.io.fload( fullfile(data_p, 'sorted_neural_data_social_gaze.mat') );
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );
ct_labels = shared_utils.io.fload( fullfile(data_p, 'celltype_labels_combined.mat') );

[spike_ts, spike_labels] = linearize_sorted( sorted );
spike_labels = apply_cell_type_labels(...
  spike_labels, ct_labels.label_mat, ct_labels.label_mat_cols );

events = add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

%%  psth

rois = {'whole_face', 'right_nonsocial_object_whole_face_matched'};

evt_mask = find( events.labels, [{'m1'}, rois] );
spk_mask = find( spike_labels, 'valid-unit' );

min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;

[psth_matrix, psth_labels, t] = compute_psth(...
  spike_ts, spike_labels, spk_mask, evts, events.labels, evt_mask, min_t, max_t, bin_width );

method = 'square_wave';
psth_matrix_smoothened = eisg.util.get_overlapping_50ms_psth_from_nonoverlapping_10ms_psth( psth_matrix, method );

%%  mean psth + peak times

n_devs = 4; % detect peak at first time above mean + std * n_devs

[mean_labs, mean_I] = retaineach( psth_labels, {'uuid', 'looks_by', 'roi'} );
psth_means = bfw.row_nanmean( psth_matrix_smoothened, mean_I );

peak_indices = find_first( above_sem(psth_means, n_devs) );
peaks = peak_indices;
peaks(peaks == 0) = nan;
peaks(~isnan(peaks)) = t(peaks(~isnan(peaks)));

%%  plot cumuluative histogram of latencies

save_p = fullfile( eisg.util.project_path, 'data/plots/latency/cum_hist', dsp3.datedir );
do_save = true;

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