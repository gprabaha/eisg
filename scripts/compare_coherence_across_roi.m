%data_p = '/Volumes/ExtSSD/social_gaze_spike_phase_coherence/';
data_p = '/gpfs/milgram/project/chang/pg496/social_gaze_sfcoherence';

roi = { 'face', 'eyes_nf', 'right_nonsocial_object' };
roi_ms = shared_utils.io.findmat( fullfile(data_p, roi) );
% roi_ms = roi_ms(100:110);

mean_each = { 'region', 'looks_by', 'roi', 'channel', 'uuid' };

[coh, coh_labels, f, t] = load_tf_measure( roi_ms, mean_each );

%%

coh_file = load( '/Volumes/ExtSSD/social_gaze_spike_phase_coherence/face/01142019_position_9.mat' );
coh_labels = fcat.from( coh_file.var );
coh = coh_file.var.coh;

%%

data_path = fullfile( eisg.util.project_path, 'processed_data');
ct_labels = load_cell_type_labels( data_path );
[uuid_I, uuids] = findall( coh_labels, 'uuid', find(coh_labels, {'valid-unit', 'maybe-valid-unit'}) );
match_I = bfw.find_combinations( ct_labels, uuids );

for i = 1:numel(uuid_I)
  if ( ~isempty(match_I{i}) )
    ct_label = cellstr( ct_labels, 'cell-type', match_I{i} );
    addsetcat( coh_labels, 'cell-type', ct_label, uuid_I{i} );
  end
end

%%

[~, transform_ind] = bfw.make_whole_face_roi( coh_labels );
coh = coh(transform_ind, :, :);

%%

plt_coh = coh;
plt_labels = coh_labels;

mask = pipe( rowmask(plt_labels) ...
  , @(m) find(plt_labels, 'm1', m) ...
  , @(m) findnone(plt_labels, '<cell-type>', m) ...
  , @(m) findnone(plt_labels, 'b', m) ...
  , @(m) findnone(plt_labels, {'face', 'eyes_nf'}, m) ...
);

plt_coh = plt_coh(mask, :, :);
plt_labels = plt_labels(mask);

pl = plotlabeled.make_spectrogram( f, t );
pl.imagesc( plt_coh, plt_labels, {'roi', 'region', 'cell-type'} );

%%

plt_coh = coh;
plt_labels = coh_labels;

mask = pipe( rowmask(plt_labels) ...
  , @(m) find(plt_labels, {'m1'}, m) ...
  , @(m) findnone(plt_labels, '<cell-type>', m) ...
  , @(m) findnone(plt_labels, 'b', m) ...
  , @(m) findnone(plt_labels, {'face', 'eyes_nf'}, m) ...
);

plt_coh = plt_coh(mask, :, :);
plt_labels = plt_labels(mask);

[coh_bands, coh_band_labels] = dsp3.get_band_means( plt_coh, plt_labels, f, {[15 25], [45 60]}, {'beta', 'gamma'});
pl = plotlabeled.make_common();
pl.x = t;
pl.lines( coh_bands, coh_band_labels, 'bands', {'roi', 'region', 'cell-type'} );

% pl.imagesc( plt_coh, plt_labels, {'roi', 'region', 'cell-type'} );