data_dir = '/Users/shadeeleazer/data/brains/intermediates';

psth_date_dirs = {'051423', '051423'};
model_date_dirs = {'051423','051423'};
data_prefixes = {'contra-shuffle-1', 'subsampled_ipsi-shuffle-1'};

psths = {};
indices_of_fixations = {};
psth_labels = {};
spike_counts = {};
mdls = {};
oob_label = '<out of bounds>';

behavior_data = load( fullfile(data_dir, 'distance_model', psth_date_dirs{1}, 'behavior.mat') );

for i = 1:numel(data_prefixes)
    fprintf( '\n %d of %d', i, numel(data_prefixes) );

    psth_prefix = data_prefixes{i};
    psth_p = fullfile( data_dir, 'distance_model', psth_date_dirs{i}, sprintf('%s-psth.mat', psth_prefix) );
    model_p = fullfile( data_dir, 'distance_model', model_date_dirs{i}, sprintf('%s-model.mat', psth_prefix ) );

    psth_data = load( psth_p );
    model_data = load( model_p );

    for ci = 1:numel(psth_data.allpsth)
        fprintf( '\n\t %d of %d', ci, numel(psth_data.allpsth) );
        cell_psth = psth_data.allpsth{ci};
        psth_index = psth_data.index_of_fixations{ci};
        self_distances = behavior_data.distances_to_m2s_eyes(psth_index);
        cell_spike_counts = psth_data.all_spike_counts_in_fixation_interval{ci};

        bin_edges = [0, 5, 10, 15, 20];
        bins = discretize( self_distances, bin_edges );
        bins(isnan(bins)) = numel( bin_edges ); %   add 1 for invalid bin

        bin_labels = string(arrayfun(@(x, y) sprintf('%d_%d', x, y), bin_edges(1:end-1), bin_edges(2:end), 'un', 0));
        bin_labels = bin_labels(:);
        bin_labels(end+1) = oob_label;

        [ic, bins] = findeachv( bins );
        mean_traces = cell2mat( cellfun(@(x) mean(cell_psth(x, :), 1), ic, 'un', 0) );
        mean_counts = cell2mat( cellfun(@(x) mean(cell_spike_counts(x)), ic, 'un', 0) );
        binned_counts = cellfun( @(x) cell_spike_counts(x), ic, 'un', 0 );

        one_psth = psth_data.spike_labels(ci, :);   
        one_psth.data_subset = categorical( data_prefixes(i) );
        psth_labels{end+1} = one_psth;

        one_psth = repmat( one_psth, numel(ic), 1 );
        one_psth.bins = bins;
        one_psth.bin_labels = bin_labels(bins);
        one_psth.traces = mean_traces;
        one_psth.counts_in_fixation_interval = mean_counts;
        one_psth.counts_in_fixation_interval_dists = binned_counts;
        psths{end+1} = one_psth;

        indices_of_fixations{end+1} = psth_data.index_of_fixations{ci};
        spike_counts{end+1} = psth_data.all_spike_counts_in_fixation_interval{ci};

    end

    mdls = [ mdls; model_data.mdls ];
end

psths = vertcat( psths{:} );
psth_labels = vertcat( psth_labels{:} );
all_spike_counts = vertcat( spike_counts{:} );

%%

[betas, ps] = cellfun( @extract_mdl_beta_coeffs, mdls, 'un', 0 );
mdl_betas = vertcat( betas{:} );
mdl_ps = vertcat( ps{:} );

 [I, beta_diff_labels] = findeach( psth_labels, {'unit_uuid'} );
 beta_diffs = nan( numel(I), 2 );
 for i = 1:numel(I)
     beta_diffs(i, :) = abs( mdl_betas(I{i}(1), :) - mdl_betas(I{i}(2), :) );
 end
 
 if ( 1 )
     low_split_thresh = 0.05;
     high_split_thresh = 0.2;
     low_diff_cells = find( beta_diffs(:, 1) < low_split_thresh );
     high_diff_cells = find( beta_diffs(:, 1) >= high_split_thresh );
 else
     bi = discretize( beta_diffs(:, 1), 2 );
     low_diff_cells = find( bi == 1 );
     high_diff_cells = find( bi == 2 );
 end

%%  mean and or median beta values across regions

clf;
%mask = psth_labels.region == 'ofc';
 mask(:) = true;
[I, id, C] = rowsets( 3, psth_labels, 'data_subset', 'region', {}, 'mask', mask );
L = plots.cellstr_join( C );
axs = plots.simplest_barsets( abs(mdl_betas(:, 1)), I, id, L ...
    , 'summary_func', @nanmedian ...
    , 'error_func', @plotlabeled.nansem ...
);
shared_utils.plot.match_ylims( axs );

%%  histograms of contra vs ipsi data

plt_x = all_spike_counts;
plt_y = behavior_data.distances_to_m2s_eyes;
%plt_y = behavior_data.m2_dist_to_m1s_eyes;
plt_y(plt_y > 20) = nan;

% [~, cell_ids] = findeach( psth_labels, 'unit_uuid' );
% target_cell_id = cell_ids.unit_uuid(5);

target_cell_id = beta_diff_labels.unit_uuid(low_diff_cells(8));

[I, C] = findeach( psth_labels, {'unit_uuid', 'data_subset', 'region', 'session', 'id_m1'} ...
    , psth_labels.unit_uuid == target_cell_id );

figure(1); clf;
axs = plots.panels( numel(I) );
for i = 1:numel(axs)
    ax = axs(i);
    cell_index = I{i};
    cell_mdl_betas = mdl_betas(cell_index, :);
    cell_mdl_ps = mdl_ps(cell_index, :);

    sig_str1 = ternary( cell_mdl_ps(1) < 0.05, '(*)', '' );
    sig_str2 = ternary( cell_mdl_ps(2) < 0.05, '(*)', '' );

    str_betas = sprintf( 'x1 = %0.3f%s, x2 = %0.3f%s' ...
        , cell_mdl_betas(1), sig_str1, cell_mdl_betas(2), sig_str2 );

    title_str = strjoin(plots.cellstr_join({C(i, :)}), ' | ');
    title_str = sprintf( '%s | %s', title_str, str_betas );
    
    index_of_fixations = indices_of_fixations{cell_index};
    subset_x = plt_x(index_of_fixations, :);
    subset_y = plt_y(index_of_fixations, :);

    scatter( ax, subset_x, subset_y, 'filled' );
    hold( ax, 'on' );

    med_x = nanmean( subset_x );
    med_y = nanmean( subset_y );
    shared_utils.plot.add_vertical_lines( ax, med_x );
    shared_utils.plot.add_horizontal_lines( ax, med_y );

    non_nan = ~isnan(subset_x) & ~isnan(subset_y);
    ps = polyfit( subset_x(non_nan), subset_y(non_nan), 1 );
    lims = get( ax, 'xlim' );
    xs = linspace( lims(1), lims(2), 10 );
    y = polyval( ps, xs );
    plot( ax, xs, y );

    title( ax, title_str );
    xlabel( ax, 'spike count' );
    ylabel( ax, 'other distance' );
end

shared_utils.plot.match_xlims( axs );
shared_utils.plot.match_ylims( axs );

%%  plot average per region

mask = psths.region == 'ofc';

[I, id, C] = rowsets( 2, psths, {'region', 'data_subset'}, 'bins', 'mask', mask, 'preserve', 2 );
[PI, PL] = plots.nest2( id, I, plots.cellstr_join(C) );
[axs, hs] = plots.simplest_linesets( psth_ts{1}, psths.traces, PI, PL ...
    , 'error_func', @(x) nan(1, size(x, 2)) ...
);

%%  plot cell summary

do_save = false;

mask = psths.region == 'bla';
mask = mask & psths.bin_labels ~= oob_label;
mask = mask & psths.unit_uuid == 'unit_uuid__1092';

%mask = psths.bin_labels ~= oob_label;

num_bins = max( psths.bins(mask) );
psth_ts = psth_data.psth_ts;

[cell_I, cell_C] = findeach( psths, {'unit_uuid', 'region'}, mask );

for i = 1:numel(cell_I)
    [pI, pC] = findeach( psths, {'region', 'data_subset', 'unit_uuid'}, cell_I{i} );
    panel_labels = arrayfun( @(x) plots.cellstr_join({pC(x, :)}), 1:size(pC, 1) );
    panel_labels = strrep( panel_labels, '_', ' ' );

    %   line plot
    %   bar plot
    num_plots_per_subset = 3;
    clf;
    axs = plots.panels( num_plots_per_subset * numel(pI) );

    term_sets = cell( numel(pI), 1 );
    for j = 1:numel(pI)
        mdl_index = psth_labels.unit_uuid == pC.unit_uuid(j) & ...
                    psth_labels.data_subset == pC.data_subset(j);
        mdl = mdls{mdl_index};

        panel_ind = pI{j};
        bI = arrayfun( @(x) intersect(panel_ind, find(psths.bins == x)), (1:num_bins)', 'un', 0 );

        line_axi = (j - 1) * num_plots_per_subset + 1;
        bar_axi = line_axi + 1;
        mdl_axi = bar_axi + 1; 

        line_ax = axs(line_axi);
        bar_ax = axs(bar_axi);
        mdl_ax = axs(mdl_axi);

        traces = cate1( cellfun(@(x) mean(psths.traces(x, :), 1), bI, 'un', 0) );
        counts = cellfun( @(x) cat(1, psths.counts_in_fixation_interval_dists{x}), bI, 'un', 0 );
        count_bin_indices = arrayfun( @(i, x) repmat(i, numel(x{1}), 1), (1:numel(counts))', counts, 'un', 0 );
        
        count_bin_x = vertcat( count_bin_indices{:} );
%         count_bin_x = count_bin_x + (rand( size(count_bin_x) ) * 2 - 1) * 0.25;
        count_bin_y = vertcat( counts{:} );

        mean_counts = cellfun( @nanmean, counts );

        terms = mdl.Coefficients.Estimate(2:3);
        r2 = mdl.Rsquared.AdjGeneralized;
        mdl_data = [ terms(:)', r2 ];
        
        mean_traces = smoothdata( traces, 2, 'smoothingfactor', 0.75 );
        h = plot( line_ax, psth_ts{1}, mean_traces );
        leg_labels = compose( "bin: %s", string(1:num_bins) );
        legend( h, leg_labels );

        bar( mdl_ax, mdl_data );
        set( mdl_ax, 'xticklabel', {'x1beta', 'x2beta', 'r2'});

        bar( bar_ax, mean_counts );
        set( bar_ax, 'xticklabel', leg_labels );
%         hold( bar_ax, 'on' );
%         scatter( bar_ax, count_bin_x, count_bin_y );
%         boxplot( bar_ax, count_bin_x, count_bin_y );

        title( bar_ax, panel_labels{j} );
        title( line_ax, panel_labels{j} );
        title( mdl_ax, panel_labels{j} );

        term_sets{j} = mdl_data;
    end

    line_axs = axs(1:num_plots_per_subset:end);
    bar_axs = axs(2:num_plots_per_subset:end);
    term_axs = axs(3:num_plots_per_subset:end);

    line_lims = cate1( arrayfun(@(x) get(x, 'ylim'), line_axs, 'un', 0) );
    bar_lims = cate1( arrayfun(@(x) get(x, 'ylim'), bar_axs, 'un', 0) );
    term_lims = cate1( arrayfun(@(x) get(x, 'ylim'), term_axs, 'un', 0) );

    set( line_axs, 'ylim', [min(line_lims(:, 1)), max(line_lims(:, 2))] );
    set( bar_axs, 'ylim', [min(bar_lims(:, 1)), max(bar_lims(:, 2))] );
    set( term_axs, 'ylim', [min(term_lims(:, 1)), max(term_lims(:, 2))] );

    split_index = term_split_table.region == cell_C.region(i);
    split_value = term_split_table.splits(split_index);

    term_sets = vertcat( term_sets{:} );
    abs_diff = abs( diff(term_sets, 1, 1) );
    is_low = abs_diff(1) < split_value;
    set_string = ternary( is_low, 'low contra ipsi diff', 'high contra ipsi diff' );

    s = get( get(axs(1), 'title'), 'string' );
    title( axs(1), sprintf('(%s) %s', set_string, s) );

    if ( do_save )
        save_p = fullfile( string(cell_C.region(i)), set_string );
        save_p = fullfile( save_p, strjoin(plots.cellstr_join({cell_C(i, :)})) );
    
       save_p = fullfile( '/Users/shadeeleazer/data/brains/plots/barpsths/contra_ipsi/GLM/' ...
           , sprintf('%s.png', save_p) );
        shared_utils.io.require_dir( fileparts(save_p) );
        saveas( gcf, save_p );
    end
end


%%  plot each cell individually


mask = psths.region == 'bla';
mask = mask & psths.bin_labels ~= oob_label;
%mask = mask & psths.unit_uuid == 'unit_uuid__1092';

[cell_I, cell_C] = findeach( psths, {'unit_uuid', 'region'}, mask );
% = cell_I(1);

psth_ts = psth_data.psth_ts;
bin_w = uniquetol( diff(psth_ts{1}) );
num_bins = numel( unique(psths.bins) );

do_save = true;
do_bar_plot = true;

for i = 1:numel(cell_I)

    [pI, pC] = findeach( psths, {'region', 'data_subset', 'unit_uuid'}, cell_I{i} );
    clf;
    axs = plots.panels( numel(pI), true );
    for j = 1:numel(pI)

        if ( do_bar_plot )
            mean_means = nan( num_bins, 1 );
        end
        for k = 1:num_bins
            bin_index = intersect( find(psths.bins == k), pI{j} );
            mean_trace = mean( psths.traces(bin_index, :), 1 );
            mean_trace = mean_trace ./ 1 / bin_w;

            if ( do_bar_plot )
                psth_t_mask = psth_ts{1} >= 0 & psth_ts{1} <= 0.25;
                mean_means(k) = mean( mean_trace(psth_t_mask) );
            else
                mean_trace = smoothdata( mean_trace, 'smoothingfactor', 0.75 );
                plot( axs(j), psth_ts{1}, mean_trace, 'displayname', strrep(bin_labels{k}, '_', ' ') );
                hold( axs(j), 'on' );
            end
        end

        pmeans = mean_means(1:end-1);
        p = polyfit( 1:numel(pmeans), pmeans, 1 );
        y = polyval( p, 1:numel(pmeans) );

        if ( do_bar_plot )
            bar( axs(j), mean_means );
            set( axs(j), 'xticklabel', strrep(bin_labels, '_', ' ') );
            xlim( axs(j), [0, numel(bin_labels)] );
            hold( axs(j), 'on' );
            plot( axs(j), 1:numel(pmeans), y );
        else
            legend;
            xlim( axs(j), [-0.1, 0.5] );
            hold( axs(j), 'on' );
        end
        
        subset_psth = psths(pI{j}(1), :);
        mdl_index = psth_labels.unit_uuid == subset_psth.unit_uuid & ...
            psth_labels.data_subset == subset_psth.data_subset;
        mdl = mdls{mdl_index};
        r2 = mdl.Rsquared.AdjGeneralized;
        curr_str = plots.cellstr_join( {pC(j, :)}, ' | ' ); 
        title( axs(j), sprintf('%s (R2 = %0.3f)', char(curr_str), r2) );
    end

    shared_utils.plot.match_ylims( axs );

    if ( ~do_bar_plot )
        shared_utils.plot.add_vertical_lines( axs, [0, 0.25] );
        ylabel( axs, 'Sp/s' );
    end
    

    if ( do_save )
       save_p = strjoin( plots.cellstr_join({cell_C(i, :)}) );
       save_p = fullfile( '/Users/shadeeleazer/data/brains/plots/barpsths/contra_ipsi/GLM/BLA/best_fit_line/norm/' ...
           , sprintf('%s.png', save_p) );
        shared_utils.io.require_dir( fileparts(save_p) );
        saveas( gcf, save_p );
    end

end

