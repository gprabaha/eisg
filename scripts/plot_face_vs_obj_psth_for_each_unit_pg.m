clear;

data_p = '/Users/prabaha/repositories/eisg/processed_data';

sorted = shared_utils.io.fload( fullfile(data_p, 'sorted_neural_data_social_gaze.mat') );
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );
ct_labels = shared_utils.io.fload( fullfile(data_p, 'celltype_labels_combined.mat') );

[spike_ts, spike_labels] = linearize_sorted( sorted );
spike_labels = apply_cell_type_labels(...
  spike_labels, ct_labels.label_mat, ct_labels.label_mat_cols );

events = add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

%% extract 10 ms nonoverlapping psth

rois = { 'whole_face', 'right_nonsocial_object_whole_face_matched' };

evt_mask = find( events.labels, [{'m1'}, rois] );
spk_mask = find( spike_labels, {'valid-unit', 'maybe-valid-unit'} );

min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;

[psth_matrix_10ms_nonoverlapping, psth_labels, t] = compute_psth(...
    spike_ts, spike_labels, spk_mask ...
  , evts, events.labels, evt_mask ...
  , min_t, max_t, bin_width );

%% plot face and object psth for each cell with smoothing for 10ms nonoverlapping psth

save_p = fullfile( eisg.util.project_path, 'plots/single_unit_psth_10ms_nonoverlapping_with_smoothing/face_vs_obj' );
do_save = true;

plt_mask = pipe( rowmask(psth_labels) ...
  , @(m) findnone(psth_labels, '<cell-type>', m) ...
);

fig_I = findall( psth_labels, {'region', 'uuid'}, plt_mask );

for i = 1:numel(fig_I)
  fi = fig_I{i};
  
  pl = plotlabeled.make_common();
  
  gcats = { 'roi' };
  pcats = {'uuid', 'cell-type', 'region'};  
  plt_labs = psth_labels(fi);
  plt_dat = psth_matrix_10ms_nonoverlapping(fi, :);
  pl.x = t;
  pl.smooth_func = @(x) eisg.util.smoothen_psth_timecourse(x);
  pl.add_smoothing = true;
  axs = pl.lines( plt_dat, plt_labs, gcats, pcats );
  
  subdirs = combs( plt_labs, {'region', 'cell-type'} );
  reg_subdir = strjoin( subdirs(1, :), '_' );
  ct_subdir = strjoin( subdirs(2, :), '_' );
  full_save_p = fullfile( save_p, reg_subdir, ct_subdir );
  
  if ( do_save )
    dsp3.req_savefig( gcf, full_save_p, prune(plt_labs), [pcats, gcats] );
  end
end

%% make 50 ms overlapping PSTH with 10 ms stride

psth_matrix_50ms_overlapping = eisg.util.get_overlapping_50ms_psth_from_nonoverlapping_10ms_psth( psth_matrix_10ms_nonoverlapping );


%% plot face and object psth for each cell without smoothing for 50ms overlapping psth

save_p = fullfile( eisg.util.project_path, 'plots/single_unit_psth_50ms_overlapping_without_smoothing/face_vs_obj' );
do_save = true;

plt_mask = pipe( rowmask(psth_labels) ...
  , @(m) findnone(psth_labels, '<cell-type>', m) ...
);

fig_I = findall( psth_labels, {'region', 'uuid'}, plt_mask );

for i = 1:numel(fig_I)
  fi = fig_I{i};
  
  pl = plotlabeled.make_common();
  
  gcats = { 'roi' };
  pcats = {'uuid', 'cell-type', 'region'};  
  plt_labs = psth_labels(fi);
  plt_dat = psth_matrix_50ms_overlapping(fi, :);
  pl.x = t;
  pl.smooth_func = @(x) eisg.util.smoothen_psth_timecourse(x);
  pl.add_smoothing = true;
  axs = pl.lines( plt_dat, plt_labs, gcats, pcats );
  
  subdirs = combs( plt_labs, {'region', 'cell-type'} );
  reg_subdir = strjoin( subdirs(1, :), '_' );
  ct_subdir = strjoin( subdirs(2, :), '_' );
  full_save_p = fullfile( save_p, reg_subdir, ct_subdir );
  
  if ( do_save )
    dsp3.req_savefig( gcf, full_save_p, prune(plt_labs), [pcats, gcats] );
  end
end

%% exract 50 ms nonoverlapping psth

min_t = -0.5;
max_t = 0.5;
bin_width = 0.05;

[psth_matrix_50ms_nonoverlapping, psth_labels, t] = compute_psth(...
    spike_ts, spike_labels, spk_mask ...
  , evts, events.labels, evt_mask ...
  , min_t, max_t, bin_width );

%% plot face and object psth for each cell without smoothing

save_p = fullfile( eisg.util.project_path, 'plots/single_unit_psth_50ms_nonoverlapping_without_smoothing/face_vs_obj' );
do_save = true;

plt_mask = pipe( rowmask(psth_labels) ...
  , @(m) findnone(psth_labels, '<cell-type>', m) ...
);

fig_I = findall( psth_labels, {'region', 'uuid'}, plt_mask );

for i = 1:numel(fig_I)
  fi = fig_I{i};
  
  pl = plotlabeled.make_common();
  
  gcats = { 'roi' };
  pcats = {'uuid', 'cell-type', 'region'};  
  plt_labs = psth_labels(fi);
  plt_dat = psth_matrix_50ms_nonoverlapping(fi, :);
  pl.x = t;
  axs = pl.lines( plt_dat, plt_labs, gcats, pcats );
  
  subdirs = combs( plt_labs, {'region', 'cell-type'} );
  reg_subdir = strjoin( subdirs(1, :), '_' );
  ct_subdir = strjoin( subdirs(2, :), '_' );
  full_save_p = fullfile( save_p, reg_subdir, ct_subdir );
  
  if ( do_save )
    dsp3.req_savefig( gcf, full_save_p, prune(plt_labs), [pcats, gcats] );
  end
end