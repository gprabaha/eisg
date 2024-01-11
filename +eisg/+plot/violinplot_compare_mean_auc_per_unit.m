function violinplot_compare_mean_auc_per_unit(...
    t, z_scored_aucs, auc_labels_f_o, region, excluded_categories,...
    unit_auc_comparison_subplots, pre_time_range, post_time_range, print_stats)

if nargin < 9
    print_stats = true;
    if nargin < 8
        pre_time_range = [-0.5 0];
        post_time_range = [0 0.5];
        if nargin < 6
            unit_auc_comparison_subplots = {'pre', 'post', 'total time'};
        end
    end
end

unit_mask = pipe( rowmask( auc_labels_f_o ), ...
    @(m) find( auc_labels_f_o, region, m ), ...
    @(m) findnone( auc_labels_f_o, excluded_categories, m )...
    );
abs_zscored_auc = abs(z_scored_aucs);

[I, C] = findall( auc_labels_f_o, {'cell-type'}, unit_mask );

num_subplots = numel(unit_auc_comparison_subplots);
for subplot_ind = 1:num_subplots
    time_period = unit_auc_comparison_subplots{subplot_ind};
    switch time_period
        case 'pre'
            binary_time_ind = t > pre_time_range(1) & t <= pre_time_range(2);
        case 'post'
            binary_time_ind = t > post_time_range(1) & t <= post_time_range(2);
        case 'total time'
            binary_time_ind = ones(size(t));
        otherwise
            error('Time period label is unknown');
    end
    plotting_struct = struct();
    for i=1:numel(C)
        plotting_struct.(C{i}) = find_unit_mean_response(abs_zscored_auc, I{i}, binary_time_ind);
    end
    subplot(1, num_subplots, subplot_ind);
    violinplot(plotting_struct, C, 'ShowData', true);
    title(time_period);
    xlabel('Cell types');
    ylabel('Mean abs z-scored AUC/unit');
    if print_stats
        [p_ranksum, ~] = ranksum(plotting_struct.narrow, plotting_struct.broad);
        [~, p_ttest] = ttest2(plotting_struct.narrow, plotting_struct.broad);
    end
    fprintf('Timeperiod: %s; Ranksum p = %0.3f; ttest2 p = %0.3f;\n', time_period, p_ranksum, p_ttest);
end
set(gcf, 'Position',  [200, 200, 900, 500]);

end


function mean_response = find_unit_mean_response(abs_zscored_auc, inds, binary_time_ind)

mean_response = mean( abs_zscored_auc(inds, binary_time_ind), 2, 'omitnan');

end