function [psth_matrix, psth_labels, t] = compute_psth(...
    spike_ts, spike_labels, spk_mask ...
  , evts, evt_labels, evt_mask ...
  , min_t, max_t, bin_width)

[spk_I, spk_C] = findall( spike_labels, {'uuid', 'session'}, spk_mask );
evt_I = bfw.find_combinations( evt_labels, spk_C(2, :), evt_mask );

[psth_matrix, t] = bfw.event_psth(...
  evts, spike_ts, evt_I, spk_I, min_t, max_t, bin_width, 'concatenate', true );
psth_labels = bfw.event_psth_labels( evt_labels, spike_labels, evt_I, spk_I, true );
assert_ispair( psth_matrix, psth_labels );

end