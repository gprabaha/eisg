function m2_pos_subset = get_m2_fixation_position(m1_t0, m1_t1, m2_t, m2_pos)

assert( numel(m2_t) == size(m2_pos, 2) );

[~, t0_ind] = min( abs(m1_t0 - m2_t) );
[~, t1_ind] = min( abs(m1_t1 - m2_t) );
m2_pos_subset = m2_pos(:, t0_ind:t1_ind);


end