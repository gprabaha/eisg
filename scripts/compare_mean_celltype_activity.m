
%%

region = 'acc';
disp(region);
% celltype_ranksum_comp(base_stats, region, base_stat_labs);
% celltype_ttest2_comp(base_stats, region, base_stat_labs);

celltype_ranksum_comp(base_stats2, region, base_stat_labs);
celltype_ttest2_comp(base_stats2, region, base_stat_labs);

region = 'bla';
disp(region);
% celltype_ranksum_comp(base_stats, region, base_stat_labs);
% celltype_ttest2_comp(base_stats, region, base_stat_labs);

celltype_ranksum_comp(base_stats2, region, base_stat_labs);
celltype_ttest2_comp(base_stats2, region, base_stat_labs);


%%

function celltype_ranksum_comp(base_stats, region, base_stat_labs)

disp('ranksum test');

[p, h] = ranksum( base_stats( find( base_stat_labs, {'b', region} ), 1), base_stats( find( base_stat_labs, {'m', region} ), 1) );
fprintf('b vs m; ranksum; p=%f, h=%d\n', p, h);

[p, h] = ranksum( base_stats( find( base_stat_labs, {'b', region} ), 1), base_stats( find( base_stat_labs, {'n', region} ), 1) );
fprintf('b vs n; ranksum; p=%f, h=%d\n', p, h);

[p, h] = ranksum( base_stats( find( base_stat_labs, {'m', region} ), 1), base_stats( find( base_stat_labs, {'n', region} ), 1) );
fprintf('m vs n; ranksum; p=%f, h=%d\n', p, h);

end

function celltype_ttest2_comp(base_stats, region, base_stat_labs)

disp('2 sample t-test test');

[h, p] = ttest2( base_stats( find( base_stat_labs, {'b', region} ), 1), base_stats( find( base_stat_labs, {'m', region} ), 1) );
fprintf('b vs m; ranksum; p=%f, h=%d\n', p, h);

[h, p] = ttest2( base_stats( find( base_stat_labs, {'b', region} ), 1), base_stats( find( base_stat_labs, {'n', region} ), 1) );
fprintf('b vs n; ranksum; p=%f, h=%d\n', p, h);

[h, p] = ttest2( base_stats( find( base_stat_labs, {'m', region} ), 1), base_stats( find( base_stat_labs, {'n', region} ), 1) );
fprintf('m vs n; ranksum; p=%f, h=%d\n', p, h);

end