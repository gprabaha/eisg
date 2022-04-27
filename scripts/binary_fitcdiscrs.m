function [fit_tables, fit_labels] = binary_fitcdiscrs(data, labels, each_I, as, bs, label_category)

validateattributes( data, {'double'}, {'column'}, mfilename, 'data' );

if ( nargin < 6 )
  label_category = [];
end

assert( numel(as) == numel(bs) );

fit_tables = cell( numel(each_I), 1 );
fit_labels = cell( size(fit_tables) );

parfor it = 1:numel(each_I)
  shared_utils.general.progress( it, numel(each_I) );
  ei = each_I{it};
  
  l = fcat;
  d = {};
  for i = 1:numel(as)
    a = as{i};
    b = bs{i};

    ind_a = find( labels, a, ei );
    ind_b = find( labels, b, ei );
    
    if ( numel(ind_a) < 4 || numel(ind_b) < 4 )
      labs = append1( fcat, labels, ei );
      p_gt = nan;
    else
      labs = append1( fcat, labels, [ind_a; ind_b] );
      p_gt = binary_fitcdiscr_perm_test( data, ind_a, ind_b );
    end
    
    d{end+1, 1} = struct( 'p', p_gt );
    maybe_set_label_category( labs, a, b, label_category );
    append( l, labs );
  end
  
  fit_tables{it} = d;
  fit_labels{it} = l;
end

fit_tables = vertcat( fit_tables{:} );
fit_labels = vertcat( fcat, fit_labels{:} );

end

function labs = maybe_set_label_category(labs, a, b, label_category)

if ( ~isempty(label_category) )
  al = strjoin( cellstr(a), ' | ' );
  bl = strjoin( cellstr(b), ' | ' );
  setcat( labs, label_category, sprintf('%s_%s', al, bl) );
end

end