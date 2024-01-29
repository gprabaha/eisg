function [dists, eye_dists, hemifield_deltas, roi_labels] = get_hemifield_center_distance(...
    pos_file, events_file, roi_file, offsets_file, all_rois, looks_by)

%   m1s_rois are "m1's regions of interest", that is, they define locations
%   on m2.
m1s_rois = roi_file.m1.rects;
m2s_rois = roi_file.m2.rects;
m2_eyes_roi = m2s_rois('eyes_nf');
m1_eyes_roi = m1s_rois('eyes_nf');

%   get center point delineating m1's hemifields
hemifield_origin = get_hemifield_zero_point( offsets_file, m2_eyes_roi );

%%

start_times = bfw.event_column( events_file, 'start_time' );

%   labels for each looking event
%   'eyes_nf' = eyes
%   'everywhere' = some other point in space not otherwise classified.
event_labels = events_file.labels;

looks_by_labels = event_labels(:, 1);
roi_labels = event_labels(:, 5);

if ( all_rois )
    roi_subset = true( size(roi_labels) );
else
    roi_subset = strcmp( roi_labels, 'eyes_nf' );
end

% roi_subset = strcmp( roi_labels, 'eyes_nf' ) | strcmp( roi_labels, 'face' );

looks_by = cellstr( looks_by );
looks_by_events = false( numel(roi_labels), 1 );

for i = 1:numel(looks_by)
    looks_by_events = looks_by_events | strcmp( looks_by_labels, looks_by{i} );
end

is_valid_event = ~isnan( start_times );

start_indices = bfw.event_column( events_file, 'start_index' );
stop_indices = bfw.event_column( events_file, 'stop_index' );

%all eye event indices
target_event_indices = find(is_valid_event & looks_by_events & roi_subset);

%%  extract mean position of each fixation event

%find all fixation start and stop times relative to eye event index
fix_event_start = start_indices(target_event_indices);
fix_event_stop = stop_indices(target_event_indices);

%cell array with all x and y coordinates for each fixation event
trial_looks_by_xy = {};
for i = 1:length(target_event_indices)
    trial_looks_by_xy{i} = pos_file.m1(:,fix_event_start(i,:):fix_event_stop(i,:));
end

%%  compute distance with respect to hemifield center

%means of all xy for m1
for i = 1:length(fix_event_start)
    allmean_xy{i} = nanmean(trial_looks_by_xy{1,i},2);
end

%delta for all means
for i = 1:size(fix_event_start)
    alldelta{i} = allmean_xy{i}' - hemifield_origin(:)';
end

hemifield_deltas = vertcat( alldelta{:} );

% preallocate an array full of nans with one element per element of `alldelta`
dists = nan( numel(alldelta), 1 );

% compute the norm of each element (`i`) of `alldelta`
for i = 1:numel(alldelta)
  dists(i) = norm( alldelta{i} );
end

if ( any(strcmp(looks_by, 'm1')) )
    assert( ~any(strcmp(looks_by, 'm2')), 'Ambiguous target roi selection' );
    target_roi = m1_eyes_roi;
else
    target_roi = m2_eyes_roi;
end

%delta for all means
for i = 1:size(fix_event_start)
    roi_center = shared_utils.rect.center( target_roi );
    alldelta{i} = allmean_xy{i}' - roi_center(:)';
end

% mean_deltas = vertcat( alldelta{:} );
% preallocate an array full of nans with one element per element of `alldelta`
eye_dists = nan( numel(alldelta), 1 );

% compute the norm of each element (`i`) of `alldelta`
for i = 1:numel(alldelta)
  eye_dists(i) = norm( alldelta{i} );
end

roi_labels = roi_labels(target_event_indices);

end