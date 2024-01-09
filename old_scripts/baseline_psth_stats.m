function [base_stats, base_stat_labs] = baseline_psth_stats(...
  evts, evt_labels, evt_mask, spike_ts, spike_labels, spk_mask, varargin)

defaults = struct();
defaults.min_t = -0.5;
defaults.max_t = 0.5;
defaults.bin_width = 0.05;
defaults.baseline_time_snippet = [];

params = shared_utils.general.parsestruct( defaults, varargin );

base_stats = [];
base_stat_labs = fcat();

[evt_I, evt_C] = findall( evt_labels, 'session', evt_mask );
for i = 1:numel(evt_I)
  fprintf( 'Session %s (%d of %d)\n', evt_C{i}, i, numel(evt_I) );
  
  sesh_evt_mask = evt_I{i};
  sesh_spk_mask = find( spike_labels, evt_C(:, i), spk_mask );
  spk_I = findall( spike_labels, 'uuid', sesh_spk_mask );
  spk_Is = shared_utils.vector.bin( spk_I, 6 );
  
  for j = 1:numel(spk_Is)
    fprintf( '\t Group (%d of %d)\n', j, numel(spk_Is) );
    si = vertcat( spk_Is{j}{:} );
    [base_psth_matrix, base_psth_labels, t] = compute_psth(...
      spike_ts, spike_labels, si, evts, evt_labels, sesh_evt_mask ...
      , params.min_t, params.max_t, params.bin_width );

    % !! IMPORTANT !!
    % !! Here time_snippet is hardcoded. Parametrize it !!
    % !! Here it is calculating the mean for 50 to 400 ms after fix

    time_snippet = t >= 0.05 & t < 0.4;
    base_psth_matrix = mean( base_psth_matrix( :, time_snippet), 2 );
    [unit_labs, unit_I] = retaineach( base_psth_labels, {'uuid', 'looks_by'} );
    sigmas = cellfun( @(x) sem_all(columnize(base_psth_matrix(x))), unit_I );
    means = cellfun( @(x) mean(columnize(base_psth_matrix(x))), unit_I );
    base_stats = [ base_stats; [means, sigmas] ];

    append( base_stat_labs, unit_labs );
  end
end

end