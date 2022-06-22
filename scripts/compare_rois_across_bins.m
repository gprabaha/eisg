clc;
clear;

data_p = fullfile( eisg.util.project_path, 'processed_data');

sorted = shared_utils.io.fload( fullfile(data_p, 'sorted_neural_data_social_gaze.mat') );
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );

%%

ct_labels = load_cell_type_labels( data_p );

[unit_spike_ts, unit_wfs, spike_labels] = linearize_sorted( sorted );
bfw.add_monk_labels( spike_labels );

[uuid_I, uuids] = findall( spike_labels, 'uuid', find(spike_labels, {'valid-unit', 'maybe-valid-unit'}) );
match_I = bfw.find_combinations( ct_labels, uuids );

for i = 1:numel(uuid_I)
  if ( ~isempty(match_I{i}) )
    ct_label = cellstr( ct_labels, 'cell-type', match_I{i} );
    addsetcat( spike_labels, 'cell-type', ct_label, uuid_I{i} );
  end
end

events = add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

%%
rois = { 'face', 'eyes_nf', 'whole_face', 'right_nonsocial_object_whole_face_matched' };

evt_mask = find( events.labels, [{'m1'}, rois] );
spk_mask = find( spike_labels, {'valid-unit', 'maybe-valid-unit'} );

min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;

[psth_matrix, psth_labels, t] = compute_psth(...
    unit_spike_ts, spike_labels, spk_mask ...
  , evts, events.labels, evt_mask ...
  , min_t, max_t, bin_width );

%% Ranksum comparison for all bins of all unit

num_bins = size( psth_matrix, 2 );

parfor i=1:num_bins
  rs_outs{i} = dsp3.ranksum( psth_matrix(:, i), psth_labels, {'uuid'}, 'whole_face', 'right_nonsocial_object_whole_face_matched' );
  % rs_outs{i} = dsp3.ranksum( psth_matrix(:, i), psth_labels, {'uuid'}, 'face', 'eyes_nf' );
  ps(:,i) = cellfun( @(x) x.p, rs_outs{i}.rs_tables );
  fprintf( 'compared bin %d of %d\n', i, num_bins );
end

%%
sig_cell_bins = zeros( size( ps ) );
sig_cell_bins(ps < 0.05) = 1;

max_consec_bins = zeros( size(sig_cell_bins, 1), 1 );
for i = 1:size(sig_cell_bins, 1)
  sig_bins = find( sig_cell_bins(i,:) );
  if ~isempty( sig_bins )
    consec_bins = 1;
    max_consec_bins(i) = 1;
    sig_bin_diff = diff( sig_bins );
    for j = 1:numel(sig_bin_diff)
      if sig_bin_diff(j) == 1
        consec_bins = consec_bins + 1;
        if consec_bins > max_consec_bins(i)
          max_consec_bins(i) = consec_bins;
        end
      else
        consec_bins = 1;
      end
    end
  end
end

