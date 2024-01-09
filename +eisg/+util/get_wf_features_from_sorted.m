function [unit_wf_features, feature_labels] = get_wf_features_from_sorted(sorted)

feature_labels = fcat();
unit_wf_features = [];

for i = 1:numel(sorted)  
  all_features = fields( sorted(i).calculated_template_features );
  
  for feature = all_features'
    unit_uuids = sorted(i).uuid;
    file_unit_features = cell2mat( sorted(i).calculated_template_features.(string( feature )) );
    
    unit_wf_features = [unit_wf_features; file_unit_features'];
    l = fcat.create( ...
        'wf_feature', feature ...
      , 'region', sorted(i).region ...
      , 'filename', sorted(i).filename ...
      , 'session', parse_session(sorted(i).filename) ...
      , 'validity', parse_validities(sorted(i).validity) ...
      , 'uuid', arrayfun(@(x) sprintf('uuid-%d', x), unit_uuids, 'un', 0) ...
      );
      append( feature_labels, l );
  end
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