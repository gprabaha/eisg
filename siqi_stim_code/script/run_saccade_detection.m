int_p = '/Volumes/external3/data/changlab/siqi/stim/intermediates';
dst_p = fullfile( int_p, 'microsaccades' );
shared_utils.io.require_dir( dst_p );

time_p = fullfile( int_p, 'aligned_raw_samples/time' );
pos_mats = shared_utils.io.findmat( fullfile(int_p, 'aligned_raw_samples/position') );
pos_fnames = shared_utils.io.filenames( pos_mats, true );

% can't parfor because `detect_microsaccades` uses a library that
% reads/writes to disk
for i = 1887:numel(pos_mats)
  
fprintf( '\n %d of %d', i, numel(pos_mats) );
  
pos_file = shared_utils.io.fload( pos_mats{i} );
t_file = shared_utils.io.fload( fullfile(time_p, pos_fnames{i}) );

first_non_nan = find( ~isnan(t_file.t), 1 );
m1_p = pos_file.m1(:, first_non_nan:end);
t = reshape( t_file.t(first_non_nan:end), [], 1 )';
fixi = do_fix_detect( m1_p, t );

fs = 1e3;
pos_deg = bfw.px2deg( m1_p );
microsaccades = detect_microsaccades( pos_deg, fs, tempdir );

isect_sample_threshold = 50;
saccade_labels = label_saccades( fixi, microsaccades, isect_sample_threshold );
fixi = make_saccade_fixation_intervals_exclusive( ...
  fixi, microsaccades(saccade_labels == 'macrosaccade', :) ...
);

is_fix = false( size(t_file.t) );
for j = 1:size(fixi, 1)
  is_fix(fixi(j, 1):fixi(j, 2)) = true;
end

save( fullfile(dst_p, pos_fnames{i}), 'microsaccades', 'saccade_labels', 'is_fix', 'fixi' );

end

%%

function fixi = make_saccade_fixation_intervals_exclusive(fixi, saccades)

need_update_intervals = true;

for i = 1:size(saccades, 1)
  if ( need_update_intervals )
    fixis = arrayfun( @(x) fixi(x, 1):fixi(x, 2), 1:size(fixi, 1), 'un', 0 );
  end
  
  sacc_interval = saccades(i, 1):saccades(i, 2);
  isects = find( cellfun(@(x) ~isempty(intersect(sacc_interval, x)), fixis) );
  
  if ( ~isempty(isects) )
    need_update_intervals = true;
  end
  
  new_fixis = [];
  for j = 1:numel(isects)
    cfixi = fixis{isects(j)};
    is_fix = true( size(cfixi) );
    [~, ia] = intersect( cfixi, sacc_interval );
    is_fix(ia) = false;
    [new_start, new_durs] = shared_utils.logical.find_islands( is_fix );
    
    for k = 1:numel(new_start)
      new_fixis(end+1, :) = [cfixi(new_start(k)), cfixi(new_start(k)+new_durs(k)-1)];
    end
  end
  
  fixi(isects, :) = [];
  fixi = [ fixi; new_fixis ];
end

[~, ord] = sort( fixi(:, 1) );
fixi = fixi(ord, :);

end

function fixi = do_fix_detect(p, t)

% repositories/eyelink
is_fix = is_fixation( p, t, 20, 10, 0.03 );
is_fix = logical( is_fix(1:numel(t)) );

[fix_starts, fix_durs] = shared_utils.logical.find_islands( is_fix );
fix_stops = fix_starts + fix_durs - 1;
fixi = [ fix_starts(:), fix_stops(:) ];

end