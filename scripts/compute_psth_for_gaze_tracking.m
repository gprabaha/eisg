function [spike_counts, ind_fix, ind_spk] = compute_psth_for_gaze_tracking(...
    spike_file, sessions, fixation_start_ts, fixation_stop_ts)

assert( numel(sessions) == numel(fixation_start_ts) );
assert( numel(fixation_start_ts) == numel(fixation_stop_ts) );

all_spike_counts_in_fixation_interval = cell( size(spike_file.spike_times) );
index_of_fixations = cell( size(all_spike_counts_in_fixation_interval) );
index_of_spikes = cell( size(all_spike_counts_in_fixation_interval) );

parfor i = 1:numel(spike_file.spike_times)

fprintf( '\n %d of %d', i, numel(spike_file.spike_times) );

cell_index = i;
         
spike_times = spike_file.spike_times{cell_index};
spike_session = char( spike_file.labels(cell_index, 'session') );

index_of_fixations_this_session = sessions == spike_session;

fix_t0s = fixation_start_ts(index_of_fixations_this_session);
fix_t1s = fixation_stop_ts(index_of_fixations_this_session);
    
all_spike_counts_in_fixation_interval{i} = count_spikes_intervals( spike_times, fix_t0s, fix_t1s );
index_of_fixations{i} = find( index_of_fixations_this_session );
index_of_spikes{i} = repmat( i, sum(index_of_fixations_this_session), 1 );

end

spike_counts = vertcat( all_spike_counts_in_fixation_interval{:} );
ind_fix = vertcat( index_of_fixations{:} );
ind_spk = vertcat( index_of_spikes{:} );

end

function c = count_spikes_in_interval(spike_times, t0, t1)

in_bounds = spike_times > t0 & spike_times <= t1;
c = sum( in_bounds );

end

function cs = count_spikes_intervals(spike_times, t0s, t1s)

cs = nan( size(t0s) );
for i = 1:numel(t0s)
  cs(i) = count_spikes_in_interval( spike_times, t0s(i), t1s(i) );
end

end