function draw_progress_bar(progress_index, total_indices, num_ticks_progress_bar)

progress_bar = repmat(' ', 1, num_ticks_progress_bar);
progress_amount = floor( (progress_index/total_indices)*num_ticks_progress_bar );
if progress_amount > 1
    progress_bar(1:progress_amount-1) = '=';
    progress_bar(progress_amount) = '>';
end
fprintf('[%s] %0.4f%% (%d/%d)\n', progress_bar, progress_index/total_indices*100, progress_index, total_indices);

end