function [fix_event_start, fix_event_stop, duration, roi_labels] = get_fixation_times(events_file, looks_by)

start_times = bfw.event_column( events_file, 'start_time' );

%   labels for each looking event
%   'eyes_nf' = eyes
%   'everywhere' = some other point in space not otherwise classified.
event_labels = events_file.labels;

looks_by_labels = event_labels(:, 1);
roi_labels = event_labels(:, 5);

looks_by = cellstr( looks_by );
looks_by_animal = false( numel(roi_labels), 1 );

for i = 1:numel(looks_by)
    looks_by_animal = looks_by_animal | strcmp( looks_by_labels, looks_by{i} );
end

is_valid_event = ~isnan( start_times );

start_times = bfw.event_column( events_file, 'start_time' );
stop_times = bfw.event_column( events_file, 'stop_time' );
duration = bfw.event_column( events_file, 'length' );

%all eye event indices
target_event_indices = find(is_valid_event & looks_by_animal);

%find all fixation start and stop times relative to eye event index
fix_event_start = start_times(target_event_indices);
fix_event_stop = stop_times(target_event_indices);
duration = duration(target_event_indices);

assert( numel(fix_event_start), size(roi_labels, 1) );

end