function valid_unit_psth = extract_valid_unit_psth(sorted_neural_data, events, varargin)

defaults = eisg.util.make_analysis_params_struct();
defaults.mask = rowmask( events.labels );

params = shared_utils.general.parsestruct( defaults, varargin );

event_data = events.events;
event_labels = events.labels;
event_key = events.event_key;

regions_per_file = { sorted_neural_data(1:end).region };
all_regions = unique( regions_per_file );

valid_unit_psth = struct();
valid_unit_psth.trial_psth = [];
valid_unit_psth.trial_psth_labels = [];
valid_unit_psth.label_categories = [];
valid_unit_psth.bin_t = [];
valid_unit_psth.bin_width = [];

categories = getcats( event_labels );
growing_psth = [];
growing_psth_labels = [];

for file_ind = 1:numel(sorted_neural_data)
    clc;
    disp('PSTH extraction progress:');
    eisg.util.draw_progress_bar(file_ind-1, numel(sorted_neural_data), params.num_ticks_in_progress_bar);
    n_units = sorted_neural_data(file_ind).n_units;
    filename = sorted_neural_data(file_ind).filename;
    dash_locs = strfind(filename, '-');
    session_id = filename(dash_locs(end)-8:dash_locs(end)-1);
    matches_session_id = find( event_labels, session_id, params.mask );
    event_start = event_data(matches_session_id, event_key('start_time'));
    labels = categorical( event_labels, getcats( event_labels ), matches_session_id );
    file_spikeinds = sorted_neural_data(file_ind).spikeindices;
    unit_number_vec = file_spikeinds(3,:);
    sampling_frequency = sorted_neural_data(file_ind).sort_params.sampling_frequency;
    min_t = -0.5;
    max_t = 0.5;
    bin_width = 0.01;
    for unit_ind = 1:n_units
        is_unit_valid = eisg.util.is_valid(sorted_neural_data, file_ind, unit_ind, params);
        if is_unit_valid
            unit_spikeinds = file_spikeinds(2, unit_number_vec == unit_ind);
            spike_times = unit_spikeinds / sampling_frequency;
            spike_times = spike_times(:);
            
            % This is for the start of roi fixation
            [psth, bin_t] = bfw.trial_psth( spike_times, event_start, min_t, max_t, bin_width );
            addtl_labels = labels;
            addtl_labels(:, end+1) = sorted_neural_data(file_ind).region;
            addtl_labels(:, end+1) = num2str(sorted_neural_data(file_ind).uuid(unit_ind));
            addtl_labels(:, end+1) = session_id;
            addtl_labels(:, end+1) = num2str(file_ind);
            addtl_labels(:, end+1) = num2str(unit_ind);
            growing_psth = [growing_psth; psth];
            growing_psth_labels = [growing_psth_labels; addtl_labels];
            
        end
    end
end
clc;
disp('PSTH extraction progress:');
eisg.util.draw_progress_bar(file_ind, numel(sorted_neural_data), params.num_ticks_in_progress_bar);

categories{end+1} = 'region';
categories{end+1} = 'uuid';
categories{end+1} = 'session_id';
categories{end+1} = 'file_ind';
categories{end+1} = 'unit_ind';

% growing_psth_labels = fcat.from( growing_psth_labels, categories );
valid_unit_psth.trial_psth = growing_psth;
valid_unit_psth.trial_psth_labels = growing_psth_labels;
valid_unit_psth.label_categories = categories;
valid_unit_psth.bin_t = bin_t;
valid_unit_psth.bin_width = bin_width;


end
