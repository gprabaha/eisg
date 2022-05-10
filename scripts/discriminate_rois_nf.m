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

%%  discriminate between pairs of rois

t_win = [0, 0.5];
t_mask = t >= t_win(1) & t < t_win(2);
t_mean = mean( psth_matrix(:, t_mask), 2 );
t_subset = psth_matrix(:, t_mask);

% discrim_x = t_mean; % use single value for each observation
discrim_x = t_subset; % use vector of values for each observation

rois = combs( psth_labels, 'roi' );
roi_pairs = bfw.pair_combination_indices( numel(rois) );

as = arrayfun( @(x) rois{x}, roi_pairs(:, 1), 'un', 0 );
bs = arrayfun( @(x) rois{x}, roi_pairs(:, 2), 'un', 0 );

% discriminate between a and b, for each a in as and b in bs
each_I = findall( psth_labels, {'uuid', 'looks_by'} );
[fit_res, fit_labels] = binary_fitcdiscrs( ...
  discrim_x, psth_labels, each_I, as, bs, 'roi' );

%%  plot results

% extract p value from classification
ps = cellfun( @(x) x.p, fit_res );
is_sig = ps < 0.05;

% proportions of significant cells per region and roi
mask = findnone( fit_labels, {'<cell-type>'} );
[prop_labs, prop_I] = retaineach(...
  fit_labels, {'region', 'looks_by', 'roi', 'cell-type'}, mask );
props = proportions( is_sig, prop_I );

pl = plotlabeled.make_common();
axs = pl.bar( props, prop_labs, {}, 'region', {'roi', 'looks_by', 'cell-type'} );
ylabel( axs(1), 'Prop Significant a vs b.' );