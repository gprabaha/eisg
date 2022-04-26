function [all_spike_ts, spike_labels] = linearize_sorted(sorted)

sr = 40e3;

spike_labels = fcat();
all_spike_ts = {};

for i = 1:numel(sorted)
  si = sorted(i).spikeindices;
  si(3, :) = sorted(i).uuid(si(3, :));  
  [unit_uuids, spike_inds] = bfw.collect_spikes( si );
  spike_ts = cellfun( @(x) (x-1)/sr, spike_inds, 'un', 0 );  
  l = fcat.create( ...
      'region', sorted(i).region ...
    , 'filename', sorted(i).filename ...
    , 'session', parse_session(sorted(i).filename) ...
    , 'validity', parse_validity(sorted(i).validity) ...
    , 'uuid', arrayfun(@(x) sprintf('uuid-%d', x), unit_uuids, 'un', 0) ...
  );
  all_spike_ts = [ all_spike_ts; spike_ts ];
  append( spike_labels, l );
end

end

function sesh = parse_session(fname)
sesh = fname(isstrprop(fname, 'digit'));
end
function str = parse_validity(valid)
str = arrayfun( @(x) ternary(x == 1, 'valid-unit', 'invalid-unit'), valid, 'un', 0 );
end