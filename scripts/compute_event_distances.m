function res = compute_event_distances(events, pos_file_list, fix_file_list, roi_file_list)

%%

pos_file_names = shared_utils.io.filenames( pos_file_list, true );
fix_file_names = shared_utils.io.filenames( fix_file_list, true );
roi_file_names = shared_utils.io.filenames( roi_file_list, true );
assert( isequal(pos_file_names(:), fix_file_names(:), roi_file_names(:)) ...
    , 'Expected matching sets of files' );

[I, C] = findeach( events.labels, 'unified_filename' );
[matched, lb] = ismember( string(C), pos_file_names );

I = I(matched);
C = C(matched);
lb = lb(matched);

pos_file_list = pos_file_list(lb);
fix_file_list = fix_file_list(lb);

m1_dist_to_m2s_eyes = nan( rows(events.events), 1 );
m2_dist_to_m1s_eyes = nan( rows(events.events), 1 );
m2_fix_props = nan( size(m2_dist_to_m1s_eyes) );

for i = 1:numel(pos_file_list)
    event_index = I{i};
    pos_file = shared_utils.io.fload( pos_file_list{i} );
    fix_file = shared_utils.io.fload( fix_file_list{i} );
    roi_file = shared_utils.io.fload( roi_file_list{i} );
    start_ind = events.events(event_index, events.event_key('start_index'));
    stop_ind = events.events(event_index, events.event_key('stop_index'));

    m2_eye = roi_file.m1.rects('eyes_nf');
    m2_eye_center = [ mean(m2_eye([1, 3])), mean(m2_eye([2, 4])) ];

    m1_eye = roi_file.m2.rects('eyes_nf');
    m1_eye_center = [ mean(m1_eye([1, 3])), mean(m1_eye([2, 4])) ];

    mean_m2_pos = arrayfun( @(s, e) nanmean(pos_file.m2(:, s:e), 2)', start_ind, stop_ind, 'un', 0 );
    mean_m2_pos = cell2mat( mean_m2_pos  );
    m2_eye_dist = vecnorm( mean_m2_pos  - m1_eye_center(:)', 2, 2 );
    m2_dist_to_m1s_eyes(event_index) = m2_eye_dist;

    mean_m1_pos = arrayfun( @(s, e) nanmean(pos_file.m1(:, s:e), 2)', start_ind, stop_ind, 'un', 0 );
    mean_m1_pos = cell2mat( mean_m1_pos );
    m1_eye_dist = vecnorm( mean_m1_pos - m2_eye_center(:)', 2, 2 );
    m1_dist_to_m2s_eyes(event_index) = m1_eye_dist;

    m2_props = arrayfun( @(s, e) nanmean(fix_file.m2(:, s:e), 2)', start_ind, stop_ind );
    m2_fix_props(event_index) = m2_props;
end

res = struct();
res.m1_dist_to_m2s_eyes = bfw.px2deg( m1_dist_to_m2s_eyes );
res.m2_dist_to_m1s_eyes = bfw.px2deg( m2_dist_to_m1s_eyes );
res.m2_fix_props = m2_fix_props;

end