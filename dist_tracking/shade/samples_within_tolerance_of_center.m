function within_tol = samples_within_tolerance_of_center(pos, dist_threshold, p_threshold)

dist = vecnorm( pos - nanmean(pos, 2), 2, 1 );
tf = dist < dist_threshold;
within_tol = pnz( tf ) >= p_threshold;

end