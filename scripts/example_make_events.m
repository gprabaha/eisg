intermediate_root = '/Volumes/ExtSSD/changlab/social_gaze_behavioral_data';
pos_dir =       fullfile( intermediate_root, 'aligned_raw_samples/position' );
time_dir =      fullfile( intermediate_root, 'aligned_raw_samples/time' );
bounds_dir =    fullfile( intermediate_root, 'aligned_raw_samples/bounds' );
fix_dir =       fullfile( intermediate_root, 'aligned_raw_samples/raw_eye_mmv_fixations' );
meta_dir =      fullfile( intermediate_root, 'meta' );
events_dir =    fullfile( intermediate_root, 'raw_events_remade' );
roi_dir =       fullfile( intermediate_root, 'rois' );

pos_mats = shared_utils.io.findmat( pos_dir );
pos_mats(is_hidden(pos_mats)) = [];

% mat_ind = 1;
mat_ind = find( contains(pos_mats, '01092019'), 1 );

pos_file = shared_utils.io.fload( pos_mats{mat_ind} );
time_file = shared_utils.io.fload( fullfile(time_dir, pos_file.unified_filename) );
meta_file = shared_utils.io.fload( fullfile(meta_dir, pos_file.unified_filename) );
events_file = shared_utils.io.fload( fullfile(events_dir, pos_file.unified_filename) );
roi_file = shared_utils.io.fload( fullfile(roi_dir, pos_file.unified_filename) );
bounds_file = shared_utils.io.fload( fullfile(bounds_dir, pos_file.unified_filename) );
fix_file = shared_utils.io.fload( fullfile(fix_dir, pos_file.unified_filename) );

%   for the rois, roi_file.m1 is the set of ROIs of *m2's* face, from
%   m1's perspective
%
%   for the bounds, bounds_file.m1 is the set of in/out of bounds indices
%   of m1's gaze position in an ROI of m2

%   find first valid time point in time file
t0 = time_file.t(find(~isnan(time_file.t), 1));

files = struct( ...
    'position', pos_file ...
    , 'time', time_file ...
    , 'rois', roi_file ...
    , 'bounds', bounds_file ...
    , 'raw_eye_mmv_fixations', fix_file ...
);

event_defaults = bfw_recording_event_defaults();
field_names = fieldnames(events_file.params);
event_params = rmfield( events_file.params, field_names(~isfield(event_defaults, field_names)) );

events = bfw.make.raw_events( files, event_defaults );

%%  computing fixations

fix_params = bfw.make.defaults.raw_fixations;
non_nan = ~isnan( time_file.t );

fix_detect = is_fixation( ...
    pos_file.m1, time_file.t(:)' ...
    , fix_params.t1, fix_params.t2, fix_params.min_duration );
fix_detect = fix_detect(1:end-1);

%%  locate cells from this session

cell_ind = find( spike_labels, meta_file.session ...
    , find(spike_labels, 'valid-unit') );

%%

function tf = is_hidden(f)

fnames = shared_utils.io.filenames( f );
tf = startsWith( fnames, '.' );

end