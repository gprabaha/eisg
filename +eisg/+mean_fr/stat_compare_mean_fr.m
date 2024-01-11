function result = stat_compare_mean_fr( psth_matrix, t, ind_a, ind_b,...
    method, time_window )

if nargin < 5
    time_window = [0 0.5];
    if nargin < 4
        method = 'ranksum';
    end
end

result = nan;
if ~(isempty(ind_a) || isempty(ind_b))
    time_vec = t > time_window(1) & t <= time_window(2);
    mean_vec_a = mean( psth_matrix(ind_a, time_vec), 2, 'omitnan' );
    mean_vec_b = mean( psth_matrix(ind_b, time_vec), 2, 'omitnan' );
    switch method
        case 'ranksum'
            [~, result] = ranksum( mean_vec_a, mean_vec_b );
        case 'ttest2'
            result = ttest2( mean_vec_a, mean_vec_b );
        otherwise
            error('!!Unfamiliar stat comparison method used in input!!');
    end
end

end