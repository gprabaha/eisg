data_dir = '/Users/shadeeleazer/Downloads/example';

run_filename = '01022019_position_1.mat';

pos_file = load( fullfile(data_dir, 'aligned_raw_samples/position', run_filename) );
time_file = load( fullfile(data_dir, 'aligned_raw_samples/time', run_filename) );
roi_file = load( fullfile(data_dir, 'rois', run_filename) );
offsets_file = load( fullfile(data_dir, 'single_origin_offsets', run_filename) );

%   contains, among other things, time-stamps of when m1 or m2 is looking
%   at their partner, and where.
events_file = load( fullfile(data_dir, 'events', run_filename) );

%   m1s_rois are "m1's regions of interest", that is, they define locations
%   on m2.
m1s_rois = roi_file.var.m1.rects;
m2s_rois = roi_file.var.m2.rects;
m1_eyes_roi = m1s_rois('eyes_nf');
m2_eyes_roi = m2s_rois('eyes_nf');
m1_r_obj_roi = m1s_rois('right_nonsocial_object');
m1_l_obj_roi = m1s_rois('left_nonsocial_object');

%   get center point delineating m1's hemifields
center_roi_m2_eyes = get_hemifield_zero_point( offsets_file.var, m2_eyes_roi );

%%

start_times = bfw.event_column( events_file.var, 'start_time' );
stop_times = bfw.event_column( events_file.var, 'stop_time' );

%   columns of event labels
event_categories = events_file.var.categories;

%   labels for each looking event
%   'eyes_nf' = eyes
%   'everywhere' = some other point in space not otherwise classified.
event_labels = events_file.var.labels;

looks_by_labels = event_labels(:, 1);
roi_labels = event_labels(:, 5);
eye_events = strcmp( roi_labels, 'eyes_nf' );
m1_events = strcmp( looks_by_labels, 'm1' );
is_valid_event = ~isnan( start_times );

start_indices = bfw.event_column( events_file.var, 'start_index' );
stop_indices = bfw.event_column( events_file.var, 'stop_index' );

first_eye_event_index = find(is_valid_event & m1_events & eye_events, 1);

%all eye event indices
eye_event_index = find(is_valid_event & m1_events & eye_events);

%%  extract mean position of each fixation event

first_fix_event_start = start_indices(first_eye_event_index);
first_fix_event_stop = stop_indices(first_eye_event_index);

%find all fixation start and stop times relative to eye event index
fix_event_start = start_indices(eye_event_index);
fix_event_stop = stop_indices(eye_event_index);


m1_xy = pos_file.var.m1(:, first_fix_event_start:first_fix_event_stop);

%cell array with all x and y coordinates for each fixation event
trial_m1_xy = {};
for i = 1:length(eye_event_index)
    trial_m1_xy{i} = pos_file.var.m1(:,fix_event_start(i,:):fix_event_stop(i,:));
end

%%  plot example rois and gaze position in roi

figure( 1 ); clf;

eye_w = m1_eyes_roi(3) - m1_eyes_roi(1);
eye_h = m1_eyes_roi(3) - m1_eyes_roi(1);

obj_w = m1_r_obj_roi(3) - m1_r_obj_roi(1);
obj_h = m1_r_obj_roi(3) - m1_r_obj_roi(1);

obj_lw = m1_l_obj_roi(3) - m1_l_obj_roi(1);
obj_lhl = m1_l_obj_roi(3) - m1_l_obj_roi(1);

rectangle( 'position', [m1_eyes_roi(1), m1_eyes_roi(2), eye_w, eye_h]  );

hold on;
rectangle( 'position', [m1_r_obj_roi(1), m1_r_obj_roi(2), obj_w, obj_h]  );
rectangle( 'position', [m1_l_obj_roi(1), m1_l_obj_roi(2), obj_lw, obj_lhl]  );

ax = gca;
set( ax, 'YDir', 'reverse' );

%   mean gaze position of eye fixation event
mean_xy = mean( m1_xy, 2 );
plot( mean_xy(1), mean_xy(2), 'k*' );

%   origin of m1's hemifield
plot( center_roi_m2_eyes(1), center_roi_m2_eyes(2), 'r*' );

%%  compute distance with respect to hemifield center

%mean_xy = mean( m1_xy, 2 );
%delta = mean_xy(:)' - center_roi_m2_eyes(:)';
%dist = norm( delta );

%means of all xy for m1
for i = 1:length(fix_event_start)
    allmean_xy{i} = mean(trial_m1_xy{1,i},2);
end

%delta for all means
for i = 1:size(fix_event_start)
    alldelta{i} = allmean_xy{i}' - center_roi_m2_eyes(:)';
end


% preallocate an array full of nans with one element per element of `alldelta`
dist = nan( numel(alldelta), 1 );

% compute the norm of each element (`i`) of `alldelta`
for i = 1:numel(alldelta)
  dist(i) = norm( alldelta{i} );
end

