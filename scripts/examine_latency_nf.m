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
bin_width = 0.05;

[psth_matrix, psth_labels, t] = compute_psth(...
  spike_ts, spike_labels, spk_mask, evts, events.labels, evt_mask, min_t, max_t, bin_width );

%%  mean psth + peak times

n_devs = 1.5; % detect peak at first time above mean + std * n_devs

[mean_labs, mean_I] = retaineach( psth_labels, {'uuid', 'looks_by', 'roi'} );
psth_means = bfw.row_nanmean( psth_matrix, mean_I );

peaks = find_first( above_std(psth_means, n_devs) );
peaks(peaks == 0) = nan;
peaks(~isnan(peaks)) = t(peaks(~isnan(peaks)));

%%  plot histogram of latencies

save_p = fullfile( eisg.util.project_path, 'data/plots/latency/hist', dsp3.datedir );
do_save = false;

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