function [base_stats, base_stat_labs] = unit_baseline_fr_across_fix(...
    psth_matrix, psth_labels, t, unit_I, baseline_time_snippet)

if nargin<6
    baseline_time_snippet = [-0.4 -0.05];
elseif baseline_time_snippet(2)<baseline_time_snippet(1)
    error('Invalid baseline time snippet!');
end

base_stats = nan(numel(unit_I), 2);
% The unit_Is are created from psth_labs to the order will be the same
base_stat_labs = retaineach(psth_labels, 'uuid');
for i = 1:numel(unit_I)
    % fprintf('Calculating for unit : %d of %d \n', i, numel(unit_I));
    if ~isempty(unit_I{i})
        time_inds = ...
            t>=baseline_time_snippet(1) & t<baseline_time_snippet(2);
        unit_mean_fr_each_trial = mean(...
            psth_matrix(unit_I{i}, time_inds), 2, 'omitnan'...
            );
        unit_mean_fr = mean(columnize(unit_mean_fr_each_trial));
        unit_fr_sem = sem_all(columnize(unit_mean_fr_each_trial));
        base_stats(i, :) = [unit_mean_fr, unit_fr_sem];
    end
end

end