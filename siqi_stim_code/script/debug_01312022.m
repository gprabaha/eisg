%%

stim_dir = 'C:\data\bfw-stim-task';
int_dir = fullfile( stim_dir, 'intermediates' );

fnames = unique( shared_utils.io.filenames(shared_utils.io.findmat(int_dir, true)) );
fnames = fnames(contains(fnames, '01312022_position_1'));
fname = sprintf( '%s.mat', fnames{1} );

roi_file = shared_utils.io.fload( fullfile(int_dir, 'rois', fname) );
pos_file = shared_utils.io.fload( fullfile(int_dir, 'aligned_raw_samples/position', fname) );
time_file = shared_utils.io.fload( fullfile(int_dir, 'aligned_raw_samples/time', fname) );
stim_file = shared_utils.io.fload( fullfile(int_dir, 'stim', fname) );

eye_roi = roi_file.m1.rects('face');
eye_roi = shared_utils.rect.expand( eye_roi, 100, 100 );

m1_pos = pos_file.m1;
in_bounds = shared_utils.rect.inside( eye_roi, m1_pos(1, :), m1_pos(2, :) );

stim_t = stim_file.stimulation_times(9);
[~, t_ind] = min( abs(stim_t - time_file.t) );

pre_s = 2e3;
post_s = 2e3;
t_range = t_ind-pre_s:t_ind+post_s;

figure(1); 
ax = gca;
cla( ax );
plot( ax, time_file.t(t_range), in_bounds(t_range) );
hold( ax, 'on' );
shared_utils.plot.add_vertical_lines( gca, time_file.t(t_ind) );