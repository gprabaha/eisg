function overlapping_psth = get_overlapping_50ms_psth_from_nonoverlapping_10ms_psth(nonoverlapping_psth)

overlapping_psth = conv2(nonoverlapping_psth', ones(5, 1)/5, 'same')';

end