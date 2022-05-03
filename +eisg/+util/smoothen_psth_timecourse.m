function smooth_psth = smoothen_psth_timecourse(timecourse)

smooth_psth = conv(timecourse, ones(10, 1)/10, 'same');

end