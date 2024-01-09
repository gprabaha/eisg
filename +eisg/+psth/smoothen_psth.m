function smooth_psth = smoothen_psth(psth)

smooth_psth = smoothdata(psth, 2, 'movmean', 10);

end