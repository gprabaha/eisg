function [all_spike_ts, all_unit_wfs, spike_labels] = linearize_sorted(sorted)

sr = 40e3;

spike_labels = fcat();
all_spike_ts = {};
all_unit_wfs = [];

for i = 1:numel(sorted)
  si = sorted(i).spikeindices;
  si(3, :) = sorted(i).uuid(si(3, :));  
  [~, spike_inds] = bfw.collect_spikes( si );
  unit_uuids = sorted(i).uuid;
  spike_ts = cellfun( @(x) (x-1)/sr, spike_inds, 'un', 0 );  
  
  unit_wfs = sorted(i).normalized_templates;
  all_unit_wfs = [ all_unit_wfs; unit_wfs ];
  
  l = fcat.create( ...
      'region', sorted(i).region ...
    , 'filename', sorted(i).filename ...
    , 'session', parse_session(sorted(i).filename) ...
    , 'validity', parse_validities(sorted(i).validity) ...
    , 'uuid', arrayfun(@(x) sprintf('uuid-%d', x), unit_uuids, 'un', 0) ...
    , 'maxchn', arrayfun( @(x) sprintf('maxchn-%d', x), sorted(i).maxchn, 'un', 0 ) ...
  );
  all_spike_ts = [ all_spike_ts; spike_ts ];
  append( spike_labels, l );
end

end

function sesh = parse_session(fname)
sesh = fname(isstrprop(fname, 'digit'));
end

function strs = parse_validities(valids)
strs = arrayfun( @parse_validity, valids, 'un', 0 );
end

function str = parse_validity(valid)
if ( isnan(valid) )
  str = 'maybe-valid-unit';
elseif ( valid == 1 )
  str = 'valid-unit';
else
  str = 'invalid-unit';
end
end