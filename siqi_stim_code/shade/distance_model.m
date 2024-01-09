%%

% 1) following the pattern for replicating subsets of labels and distances,
%   create a new cell array for all delta positions for each psth; soething like
%   `all_psth_delta_pos`
% 2) split rows of psth_by_region{i} depending on whether the x coordinates
%       of the corresponding elements of `all_psth_delta_pos` are positive
%       or negative.
%       pos_psth_elements = all_psth_delta_pos{...} ...
%       neg_psth_elements = all_psth_delta_pos{...} ...
%    the code should be similar to how we split fixations by close / med
%    distance

%%

data_dir = '/Users/shadeeleazer/data/brains/intermediates';

pos_files = shared_utils.io.findmat( fullfile(data_dir, 'aligned_raw_samples/position') );

matched_files = bfw.matched_files( pos_files ...
    , fullfile(data_dir, 'aligned_raw_samples/time') ...
    , fullfile(data_dir, 'remade_events') ...
    , fullfile(data_dir, 'rois') ...
    , fullfile(data_dir, 'single_origin_offsets') ... 
    , fullfile(data_dir, 'meta'));
%%

spike_file = shared_utils.io.fload( fullfile(data_dir, 'relabeled_cells.mat') );
bfw.add_monk_labels( spike_file.labels );

spike_labels = fcat.totable( spike_file.labels );

%%
% create loop to create matrix for distances for each run

all_m1_labels = cell( size(matched_files, 1), 1 );
all_fixation_start_ts = cell( size(all_m1_labels) );
all_fixation_stop_ts = cell( size(all_m1_labels) );
all_distances_to_self_origin = cell( size(all_m1_labels) );
all_distances_to_m2s_eyes = cell( size(all_m1_labels) );
all_m2_distances_to_m1s_eyes = cell( size(all_m1_labels) );
all_m1_hemifield_origin_deltas = cell( size(all_m1_labels) );
all_m2_valid = cell( size(all_m1_labels) ); 
all_m1_durations = cell( size(all_m1_labels) );

%    Make another aggregate array for the fixation time stamps and
%   append to it, like we did for the distances

%includes elements for each run
parfor i = 1:length(matched_files)
    %%

    fprintf( '\nProcessing: (%d of %d)', i, length(matched_files));

    one_file = matched_files(i,:);

    pos_file = shared_utils.io.fload( one_file{1} );
    time_file = shared_utils.io.fload( one_file{2} );
    events_file = shared_utils.io.fload( one_file{3} );
    roi_file = shared_utils.io.fload( one_file{4} );
    offsets_file = shared_utils.io.fload( one_file{5} );
    meta_file = shared_utils.io.fload( one_file{6} );

    %%
    
    [m1_fix_ts, m1_stop_ts, m1_durations] = get_fixation_times( events_file, {'m1', 'mutual'} );
    [dist_to_self_origin, dist_to_m2s_eyes, m1_hemifield_origin_delta_pos, m1_roi_labels] = ...
        get_hemifield_center_distance( pos_file, events_file, roi_file, offsets_file, true, {'m1', 'mutual'} );

    m2_dist = nan( size(m1_fix_ts) );
    m2_would_pass = false( size(m1_fix_ts) );
    for j = 1:numel(m1_fix_ts)
        [m2_dist(j), m2_would_pass(j)] = m2_distance_during_m1_fixation( ...
            m1_fix_ts(j), m1_stop_ts(j), time_file, pos_file, roi_file );
    end


    %%

    %   append to arrays
    all_m1_labels{i} = repmat( struct2cell( meta_file )', numel(m1_fix_ts), 1 );
    assert( iscell(all_m1_labels{i}) );

    all_fixation_start_ts(i) = { m1_fix_ts };
    all_fixation_stop_ts(i) = { m1_stop_ts };
    all_m1_durations{i} = m1_durations;
    all_distances_to_self_origin{i} = dist_to_self_origin;
    all_distances_to_m2s_eyes{i} = dist_to_m2s_eyes;
    all_m1_hemifield_origin_deltas{i} = m1_hemifield_origin_delta_pos;
    all_m2_distances_to_m1s_eyes{i} = m2_dist;
    all_m2_valid{i} = m2_would_pass;
end
%% 

%%  create matrix of datas and labels by vertically concenating sub-matrices

mat_labs = vertcat( all_m1_labels{:} );
[id_m1, id_m2] = extract_m1_m2_ids( mat_labs(:, 3) );

table_labs = array2table( ...
  [mat_labs, id_m1], 'VariableNames' ...
  , {'unified_filename', 'date', 'session', 'run', 'task_type', 'run_number', 'id_m1'} );

fixation_start_ts = vertcat( all_fixation_start_ts{:} );
fixation_stop_ts = vertcat( all_fixation_stop_ts{:} );
distances_to_self_origin = vertcat( all_distances_to_self_origin{:} );
distances_to_m2s_eyes = bfw.px2deg( vertcat(all_distances_to_m2s_eyes{:}) );
m2_dist_to_m1s_eyes = vertcat( all_m2_distances_to_m1s_eyes{:} );
m1_hemifield_origin_delta_pos = bfw.px2deg( vertcat( all_m1_hemifield_origin_deltas{:} ) );
m2_dist_is_valid = vertcat( all_m2_valid{:} );
m1_fix_durations = vertcat( all_m1_durations{:} );

[I, C] = findeach( id_m1, 1 );
for i = 1:numel(I)
   if ( strcmp(C{i}, 'm1_lynch') )
        %   lynch his chamber's on the left. 
      m1_hemifield_origin_delta_pos(I{i}) = m1_hemifield_origin_delta_pos(I{i});

    elseif ( strcmp(C{i}, 'm1_kuro') )
        %Kuro's chamber is on the right
        %kuro's sign needed to be flipped bc their chambers are on opposite
        %sides. So kuro's contra is left while lynch's is right
       m1_hemifield_origin_delta_pos(I{i}) = -m1_hemifield_origin_delta_pos(I{i});
  else
       error( 'Expected lynch or kuro.' );
   end
end

%%  make trial table with inputs for distance model

sessions = categorical( table_labs.session );

trial_table = table( ...
    distances_to_m2s_eyes ...
    , m2_dist_to_m1s_eyes ...
    , m2_dist_is_valid ...
    , m1_fix_durations ...
    , fixation_start_ts ...
    , fixation_stop_ts ...
    , m1_hemifield_origin_delta_pos ...
    , sessions );

%%

if ( 1 )
    save_p = fullfile( data_dir, 'distance_model', dsp3.datedir, 'behavior.mat' );
    shared_utils.io.require_dir( fileparts(save_p) );
    save( save_p, '-v7.3' ...
        , 'fixation_start_ts', 'fixation_stop_ts', 'distances_to_self_origin', 'distances_to_m2s_eyes' ...
        , 'm2_dist_to_m1s_eyes', 'm1_hemifield_origin_delta_pos', 'm2_dist_is_valid' ...
        , 'm1_fix_durations', 'table_labs' );
end

%%  plot psth of one cell at a time

cell_index = 2;         
spike_times = spike_file.spike_times{cell_index};
spike_session = spike_labels.session(cell_index);

index_of_fixations_this_session = categorical( table_labs.session) == spike_session;
fix_t0s = fixation_start_ts(index_of_fixations_this_session);

min_t = -0.5;
max_t = 0.5;
bin_width = 0.05;
[psth, bin_t] = bfw.trial_psth( spike_times(:), fix_t0s(:), min_t, max_t, bin_width );

figure( 1 );
clf;
% plot sum of spike counts within each time bin, over all trials
plot( bin_t, sum(psth, 1) );

%%  compute shuffled partitions of contra / ipsi

I = findeach( table_labs, 'session' );

num_iters = 1000;
contra_ipsi_labels = cell( num_iters, 1 );

for it = 1:num_iters
    contra_ipsi_label = nan( size(table_labs, 1), 1 );
    for i = 1:numel(I)
        num_contra_ipsi = floor( numel(I{i}) * 0.5 );
        all_indices = randsample( numel(I{i}), num_contra_ipsi * 2 );
        contra_inds = I{i}(all_indices(1:num_contra_ipsi));
        ipsi_inds = I{i}(all_indices(num_contra_ipsi+1:end));
        %   arbitrarily assign contra as 1 and ipsi as 2
        contra_ipsi_label(contra_inds) = 1;
        contra_ipsi_label(ipsi_inds) = 2;
    end
    
    contra_ipsi_labels{it} = contra_ipsi_label;
end

%%  load shuffled contra ipsi labels if not already in memory

load( fullfile(data_dir, 'distance_model', 'shuffled_contra_ipsi_indices', 'subsets.mat') );

%% compute psth for all cells

sessions = trial_table.sessions;
mdl_prefix = 'subsampled_ipsi';
% mdl_prefix = 'contra';

do_shuffle = true;

if ( do_shuffle )
    num_iters = 100;
else
    num_iters = 1;
end

for mdl_iteration = 1:1

%%

if ( do_shuffle )
    is_contra = contains( mdl_prefix, 'contra' );
    if ( is_contra )
        search_label = 1;
    else
        search_label = 2;
    end
    contra_ipsi_mask = contra_ipsi_labels{mdl_iteration} == search_label;
    shuffle_prefix = sprintf( '-shuffle-%d', mdl_iteration );
else
    shuffle_prefix = '';
    contra_ipsi_mask = make_real_contra_ipsi_mask( m1_hemifield_origin_delta_pos, sessions, mdl_prefix );
end

%%  compute psth using some subset of trials given by contra_ipsi_mask

include_psth_timecourse = false;

[allpsth, psth_ts, all_spike_counts_in_fixation_interval, index_of_fixations] = ...
    compute_contra_ipsi_psth( spike_file, sessions ...
    , contra_ipsi_mask, fixation_start_ts, fixation_stop_ts, include_psth_timecourse );

assert( numel(allpsth) == numel(spike_file.spike_times));

if ( 1 )
    save_p = fullfile( data_dir, 'distance_model', dsp3.datedir, sprintf('%s-psth.mat', mdl_prefix) );
    shared_utils.io.require_dir( fileparts(save_p) );
    save( save_p, '-v7.3', "allpsth", "psth_ts", "all_spike_counts_in_fixation_interval", "spike_labels", "index_of_fixations", "table_labs" );
end

%%  run distance model for all cells

%   for each cell, compute the distance model and obtain the model output
%   `mdl`.

model_table = trial_table;
model_table.contra_ipsi_mask = contra_ipsi_mask;

use_discretized = false;

discretized_prefix = '';
if ( use_discretized )
    discretized_prefix = 'disc-';
end

[mdls, cell_labels, fit_function_name] = run_distance_model( spike_file, model_table, all_spike_counts_in_fixation_interval, use_discretized );

if ( 1 )
    file_name = sprintf( '%s%s%s-model.mat', discretized_prefix, mdl_prefix, shuffle_prefix );
    save_p = fullfile( data_dir, 'distance_model', dsp3.datedir, file_name );
    shared_utils.io.require_dir( fileparts(save_p) );
    save( save_p, "mdls", "cell_labels", "spike_labels", "fit_function_name" );
end

end

%% histogram to see ipsi or contra fixation bias 
[I, C] = findeach( id_m1, 1 );
for i = 1:numel(I)
    ax = subplot( 2, 1, i );

   % edges = -20:1/2:20;
%    
    
%    
    histogram(m1_hemifield_origin_delta_pos(I{i}),edges);
    [~, ~, bin_indices] = histcounts(m1_hemifield_origin_delta_pos(I{i}),edges);

    hold on;

   % ylim( [0, 500] );

    med = nanmedian(m1_hemifield_origin_delta_pos(I{i}));
    shared_utils.plot.add_vertical_lines( gca, med );

    is_significant = signrank(m1_hemifield_origin_delta_pos(I{i})) < 0.05;

     plot_str = sprintf('Median = %0.3f', med);
    if ( is_significant )
        plot_str = sprintf( '%s (*)', plot_str );
    end

    text( ax, med, 490, plot_str );
    
    title( ax, strrep(C{i}, '_', ' ') );
    xlabel('degrees')
    ylabel('frequency')
end 

%%  plot frequencies of cell type labels separately for each region
%label 1-self cell, label2-other, label3-both, label4-neither cell

figure(1);
clf;

[I, C] = findeach( spike_labels, 'region' );
for i = 1:numel(I)

    ax = subplot(4,1,i);
    cell_labels_this_region = cell_labels(I{i});
label1 = sum(cell_labels_this_region ==1);
label2 = sum(cell_labels_this_region ==2);
label3 = sum(cell_labels_this_region ==3);
label4 = sum(cell_labels_this_region ==4);
 
 bar([label1, label2, label3, label4]);
 ylim([0 225]);
 title( C.region(i) );
end

%%  plot pie chart for frequencies of cell types

mdl_prefix = 'all_contra_ipsi';
date_dir = '040223';
model_config_p = fullfile( data_dir, 'distance_model', date_dir, sprintf('%s-model.mat', mdl_prefix) );
model_outs = load( model_config_p );

figure(3);

[I, C] = findeach( model_outs.spike_labels, 'region' );
for i = 1:numel(I)
    ax = subplot(4,1,i);
    cell_labels_this_region = model_outs.cell_labels(I{i});
    label1 = sum(cell_labels_this_region ==1);
    label2 = sum(cell_labels_this_region ==2);
    label3 = sum(cell_labels_this_region ==3);
    label4 = sum(cell_labels_this_region ==4);

pie([label1, label2, label3, label4]);
title( 'All fixations', C.region(i) );
labels = {'Self','Other', 'Both', 'Neither'};
lgd = legend(labels);

end

%%  compare frequencies of significant cells for contra vs ipsi

mdl_prefix1 = 'contra';
mdl_prefix2 = 'subsampled_ipsi';

date_dir = '041423';
model_outs1 = load( ...
    fullfile( data_dir, 'distance_model', date_dir, sprintf('%s-model.mat', mdl_prefix1)));
model_outs2 = load( ...
    fullfile( data_dir, 'distance_model', date_dir, sprintf('%s-model.mat', mdl_prefix2)));

[I, C] = findeach( model_outs1.spike_labels, 'region' );
p_significant_vs_nonsignificant = nan( numel(I), 1 );
p_factor_proportions = nan( numel(I), 1 );
for i = 1:numel(I)
    cell_labels_this_region1 = model_outs1.cell_labels(I{i});

    label1 = sum(cell_labels_this_region1 ==1);
    label2 = sum(cell_labels_this_region1 ==2);
    label3 = sum(cell_labels_this_region1 ==3);
    label4 = sum(cell_labels_this_region1 ==4);

    num_significant1 = sum(label1 + label2 + label3);
    total_num1 = sum(label1 + label2 + label3 + label4);

    cell_labels_this_region1 = model_outs2.cell_labels(I{i});

    model_outs2_label1 = sum(cell_labels_this_region1 ==1);
    model_outs2_label2 = sum(cell_labels_this_region1 ==2);
    model_outs2_label3 = sum(cell_labels_this_region1 ==3);
    model_outs2_label4 = sum(cell_labels_this_region1 ==4);

    num_significant2 = sum(model_outs2_label1 + model_outs2_label2 + model_outs2_label3);

    freq_labels = fcat.create( 'factor_type', {'self', 'other', 'both'} );
    addcat( freq_labels, 'model_type' )
    repset( freq_labels, 'model_type', {mdl_prefix1});
    repset( freq_labels, 'model_type', {mdl_prefix1, mdl_prefix2} );
    freqs = [label1; label2; label3];
    freqs = [label1; label2; label3; model_outs2_label1; model_outs2_label2; model_outs2_label3];
    chi2_info = dsp3.chi2_tabular_frequencies( freqs, freq_labels, {}, 'model_type', 'factor_type' );
    

    X1 = num_significant1;
    X2 = num_significant2; 
    X = [X1, X2];
    N = [total_num1, total_num1];

    [~,p, chi2stat] = prop_test(X , N, false);
    p_significant_vs_nonsignificant(i) = p;
    p_factor_proportions(i) = chi2_info.p;
end

C.p_significant_vs_nonsignificant = p_significant_vs_nonsignificant;
C.p_factor_proportions = p_factor_proportions;



    %   X(1) should be the total number of signifant cells, summed over
    %   distance model factor labels, for the first model configuration
    %   X(2) should be the same for the other model configuration
    %   N(1:2) should be the total number of cells, which will be the same
    %   between regions

%% same as above just get the r^2 adjusted generalized values 

mdl_names = { 'contra', 'subsampled_ipsi'};
date_dirs = {'051423', '051423'}; %datedirs for shuffled data
date_dirs = {'040723', '041423'};

shuffled = false;
shuffled_n = 100;
if ( shuffled )
    mdl_names_contra = arrayfun( @(x) string(sprintf('contra-shuffle-%d', x)), 1:shuffled_n );
    mdl_names_ipsi = arrayfun( @(x) string(sprintf('subsampled_ipsi-shuffle-%d', x)), 1:shuffled_n );
    mdl_names = string( mdl_names );
    date_dirs = string( date_dirs );
    mdl_names = [ mdl_names, mdl_names_contra, mdl_names_ipsi ];
    date_dirs = [ date_dirs, repmat("051323", 1, numel(mdl_names_contra) + numel(mdl_names_ipsi)) ];
    assert( numel(mdl_names) == numel(date_dirs) );
end

all_radj_gen = {};
all_betas = {};
all_ps = {};
all_mdl_labels = {};

for idx = 1:numel(mdl_names)
    fprintf( '\n %d of %d', idx, numel(mdl_names) );
    mdl_name = mdl_names{idx};

    output_mdls = load( ...
        fullfile( data_dir, 'distance_model', date_dirs{idx}, sprintf('%s-model.mat', mdl_name)));
    mdl_labels = output_mdls.spike_labels;
    mdl_labels.data_subset = repmat( categorical({mdl_name}), size(mdl_labels, 1), 1 );
    all_mdl_labels{end+1} = mdl_labels;
    
    for i = 1: numel(output_mdls.mdls)
        radj_gen = output_mdls.mdls{i}.Rsquared.AdjGeneralized;
        all_radj_gen{end+1} = radj_gen;
    
        [betas, ps] = extract_mdl_beta_coeffs( output_mdls.mdls{i} );
        all_betas{end+1} = betas;
        all_ps{end+1} = ps;
    end
end

all_mdl_labels = vertcat( all_mdl_labels{:} );
all_betas = vertcat( all_betas{:} );
all_ps = vertcat( all_ps{:} );

all_radj_gen = cell2mat(all_radj_gen)';
mean_all_contra_radj_gen = mean(all_radj_gen);

%%  look at sets of cells with positive or negative betas for contra or ipsi

desired_sign = 1;
% desired_subset = 'contra';
desired_subset = 'subsampled_ispi';
term_index = 1;

signs = sign( all_betas(:, term_index) );
mask = all_mdl_labels.data_subset == desired_subset & signs == desired_sign;
uuids = all_mdl_labels(mask, {'region', 'unit_uuid'});

%%  distribution of pairs of signs of model terms between contra and ipsi
%of the different sign pairs should also do the |contra| - |ipsi| and 
% |contra| / |ipsi| (could be better for normalizing)

pairs = { {'contra', 'subsampled_ipsi'} };
signs = { [1, 1], [-1, 1], [-1, -1], [1, -1] };
term_index = 1;


[I, C] = findeach( all_mdl_labels, 'region' );
rows = {};

for idx = 1:numel(I)
for i = 1:numel(pairs)
    for j = 1:numel(signs)
        ind_a = all_mdl_labels.data_subset == pairs{i}{1};
        ind_b = all_mdl_labels.data_subset == pairs{i}{2};

        ind_a = intersect( I{idx}, find(ind_a) );
        ind_b = intersect( I{idx}, find(ind_b) );

        matches_a = sign( all_betas(ind_a, term_index) ) == signs{j}(1);
        matches_b = sign( all_betas(ind_b, term_index) ) == signs{j}(2);
        num_matches = sum( matches_a & matches_b ) / numel( ind_a );

        row = table( num_matches, strjoin(string(pairs{i}), ' and '), strjoin(string(signs{j}), '/') ...
            , 'va', {'num_matches', 'pairs', 'signs'} );
        row = [ row, C(idx, :) ];
        rows{end+1, 1} = row;
    end
end
end



rows = vertcat( rows{:} );

[I, id, C] = rowsets( 3, rows, {'region'}, {'signs'}, {'pairs'} );
L = plots.cellstr_join( C );
[axs, hs, xs] = plots.simplest_barsets( rows.num_matches, I, id, L );


%%  compute beta differences for shuffled and real contra/ipsi subsets

pairs = { {'contra', 'subsampled_ipsi'} };
term_index = 1;
betas = all_betas(:, term_index);
abs_betas = abs( betas );


[I, C] = findeach( all_mdl_labels, {'region', 'unit_uuid'} );
rows = {};

for idx = 1:numel(I)
    fprintf( '\n %d of %d', idx, numel(I) );
    for i = 1:numel(pairs)
        for j = 1:shuffled_n+1
            label_a = pairs{i}{1};
            label_b = pairs{i}{2};
            if ( j <= shuffled_n )
                is_null = true;
                label_a = sprintf( '%s-shuffle-%d', label_a, j );
                label_b = sprintf( '%s-shuffle-%d', label_b, j );
            else
                is_null = false;
            end

            ind_a = all_mdl_labels.data_subset == label_a;
            ind_b = all_mdl_labels.data_subset == label_b;

            ind_a = intersect( I{idx}, find(ind_a) );
            ind_b = intersect( I{idx}, find(ind_b) );

            diffs = abs_betas(ind_a) - abs_betas(ind_b);
            sign_a = sign( betas(ind_a) );
            sign_b = sign( betas(ind_b) );

            row = table( diffs, sign_a, sign_b, is_null, 'va', {'beta_difference', 'sign_a', 'sign_b', 'null'} );
            row = [row, C(idx, :)];
            rows{end+1, 1} = row;
        end
    end
end

rows = vertcat( rows{:} );

%%  plot proportions of significant contra minus ipsi beta values, across cells

[I, C] = findeach( rows, {'unit_uuid', 'region'} );
ps = nan( size(C, 1), 1 );
real_diffs = nan( size(ps) );
sign_labels = cell( size(ps) );

ax = gca; 
cla( ax );
for i = 1:numel(I)
    is_null = rows.null(I{i});
    null_indices = I{i}(is_null);
    real_index = I{i}(~is_null);

    null_betas = rows.beta_difference(null_indices);
    real_beta = rows.beta_difference(real_index);

    if ( real_beta >= 0 )
        num_null_gt_real = sum( null_betas > real_beta );
    else
        num_null_gt_real = sum( null_betas < real_beta );
    end

    p_null_gt_real = num_null_gt_real / sum( is_null );
    ps(i) = p_null_gt_real;
    real_diffs(i) = real_beta;
    sign_labels{i} = sprintf( '%d_%d', rows.sign_a(real_index), rows.sign_b(real_index) );
end

is_sig = double( ps < 0.05 );
cell_labels = C;
cell_labels.is_sig = is_sig;
cell_labels.beta_diff = real_diffs;
cell_labels.sign_labels = sign_labels;

[I, id, C] = rowsets( 3, C, {}, {'region'}, {});
axs = plots.simplest_barsets( is_sig, I, id, plots.cellstr_join(C) ...
    , 'error_func', @(x) nan ...
);

%% Real beta difference for significant cells

beta_diff_mask = cell_labels.is_sig == 1;
% beta_diff_mask(:) = true;
beta_diff_mask = beta_diff_mask & ...
     strcmp( cell_labels.sign_labels, '-1_-1' );

[I, id, C] = rowsets( 3, cell_labels, {}, {'region'}, {} ...
    , 'mask', beta_diff_mask ...
);

clf;
axs = plots.simplest_barsets( cell_labels.beta_diff, I, id, plots.cellstr_join(C) );
ylabel( axs(1), 'Real beta difference for significant cells' );

%% of the different sign pairs should also do the |contra| - |ipsi| and 
% |contra| / |ipsi| (could be better for normalizing)


term_index = 1; 
contra_ind = all_mdl_labels.data_subset == 'contra';
ipsi_ind = all_mdl_labels.data_subset == 'subsampled_ipsi';

one_set_labels = all_mdl_labels(contra_ind, :);

%subset_ind = contra_ind & ipsi_ind & region_ind;

%[I, C] = findeach( all_mdl_labels, {'region', 'data_subset'}, subset_ind );


%when both betas are negative
   contra_betas = all_betas(contra_ind, term_index);
   negative_contra = contra_betas <0;
   positive_contra = contra_betas > 0;
   ipsi_betas = all_betas(ipsi_ind, term_index);
   negative_ipsi = ipsi_betas <0;
   positive_ipsi = ipsi_betas >0; 
   both_negative = negative_ipsi & negative_contra ==1; 
   contra_negative_betas = contra_betas( both_negative);
   ipsi_negative_betas = ipsi_betas(both_negative);
   ipsi_pos_betas = ipsi_betas(positive_ipsi);
   only_contra_negative = negative_contra ==1 & positive_ipsi==1;
   ipsi_neg_contra_pos = negative_ipsi ==1 & positive_contra ==1; 
   contra_pos_ipsi_pos = positive_contra ==1 & positive_ipsi ==1;




        positive_ipsi_unit_ids = one_set_labels.unit_uuid(positive_ipsi);
        negative_ipsi_unit_ids = one_set_labels.unit_uuid(negative_ipsi);
        neg_neg_unit_ids = one_set_labels.unit_uuid(both_negative);
        neg_pos_unit_ids = one_set_labels.unit_uuid(only_contra_negative); 
        pos_neg_unit_ids = one_set_labels.unit_uuid(ipsi_neg_contra_pos);
        pos_pos_unit_ids = one_set_labels.unit_uuid(contra_pos_ipsi_pos);

        



%% relationships between contra / ipsi beta values (- or /)

assert( isequal(cell_labels.unit_uuid, one_set_labels.unit_uuid) );

op = 'abs_diff';
%op = 'abs_proportion';
%summary_func = @mean;
summary_func = @median;
err_func = @plotlabeled.nansem;
use_only_sig = true;
do_bar = false;

switch ( op )
    case 'abs_diff'
        abscontra = abs(contra_betas);
       absipsi = abs(ipsi_betas);
       contra_ipsi_metric = abscontra - absipsi;
    case 'abs_proportion'
         abscontra = abs(contra_betas);
       absipsi = abs(ipsi_betas);
        contra_ipsi_metric = abscontra ./ absipsi;
    otherwise
        error( 'Unrecognized op "%s".', op );
end

sign_labels = cell( size(one_set_labels, 1), 1 );
sign_labels(both_negative) = {'-/-'};
sign_labels(positive_contra & negative_ipsi) = {'+/-'};
sign_labels(positive_contra & positive_ipsi) = {'+/+'};
sign_labels(negative_contra & positive_ipsi) = {'-/+'};
assert( ~any(cellfun('isempty', sign_labels)) );

one_set_labels.sign = sign_labels;
one_set_labels.is_sig = cell_labels.is_sig;

if ( use_only_sig )
    mask = find( one_set_labels.is_sig );
else
    mask = rowmask( one_set_labels );
end

[I, reg_C] = findeach( one_set_labels, {'region'}, mask );
figure(2); clf;
axs = plots.panels( numel(I) );
reg_C.signrank_ps = cell( numel(I), 1 );

for i = 1:numel(I)
    [sign_I, sign_C] = findeach( one_set_labels, {'sign', 'region'}, I{i} );
    [~, ord] = sort( sign_C.sign );
    sign_I = sign_I(ord);
    sign_C = sign_C(ord, :);

    raw_metrics = cellfun( @(x) contra_ipsi_metric(x), sign_I, 'un', 0 );
    signrank_ps = cellfun( @signrank, raw_metrics );

    if ( do_bar )
        mean_metric = cellfun( @(x) summary_func(contra_ipsi_metric(x)), sign_I );
        err_metric = cellfun( @(x) err_func(contra_ipsi_metric(x)), sign_I );
    
        bar( axs(i), mean_metric );
        set( axs(i), 'xticklabel', sign_C.sign );
    
        sig_signs = find( signrank_ps < 0.05 );
        hold( axs(i), 'on' );
        scatter( axs(i), sig_signs, mean_metric(sig_signs) + mean_metric(sig_signs) * 0.05 );
         %plots.barerrs( axs(i), (1:numel(err_metric))', mean_metric, err_metric );
    else
        [group, v] = ungroupi( sign_I );
        X = contra_ipsi_metric(v);
        subplot( axs(i) );
        violinplot( X, group );

        sig_signs = find( signrank_ps < 0.05 );
        sign_labs = sign_C.sign;
        sign_labs(sig_signs) = compose( '%s (*)', string(sign_labs(sig_signs)) );
        set( axs(i), 'xticklabel', sign_labs );
    end

    title( axs(i), reg_C.region(i) );
    reg_C.signrank_ps{i} = signrank_ps;
end

if ( do_bar )
    shared_utils.plot.match_ylims( axs );
else
    shared_utils.plot.set_ylims( axs, [-0.1, 0.1] );
end         

[I, sign_C] = findeach( one_set_labels, {'sign'}, mask );
sign_C.kw_p = nan( numel(I), 1 );

for i = 1:numel(I)
    reg_I = findeach( one_set_labels, 'region', I{i} );
    [group, v] = ungroupi( reg_I );
    X = contra_ipsi_metric(v);
    [p, stats] = kruskalwallis( X, group, 'off' );
    sign_C.kw_p(i) = p;
end

sig_ps = find( sign_C.kw_p < 0.05 );
fprintf( 'Signs that are significant across regions\n:');
curr_labs = get( axs(1), 'xticklabels' );
disp( curr_labs(sig_ps) );



% for i = 1:numel(axs)
%     curr_labs = get( axs(i), 'xticklabels' );
%     curr_labs(sig_ps) = compose( '%s (* across regions)', string(curr_labs(sig_ps)) );
%     set( axs(i), 'xticklabels', curr_labs );
% end

%% paired t-test between contra and ipsi r^2 adj generalized values 

x_pairedt = all_contra_radj_gen;
y_pairedt = all_subsampled_ipsi_radj_cell_gen';

[h,p] = ttest(x_pairedt,y_pairedt);

%h = 1 and p = 0.0033

y_all_ipsi = all_ipsi_radj_gen;

[h1,p1] = ttest(x_pairedt, y_all_ipsi); 
%h1 = 1 and p = 0.0134


%% histogram of frequency of r^2 for contra and ipsi
figure(2)

subplot(3,1,1)
histogram(all_contra_ipsi_radj_gen)

subplot(3,1,2)
histogram(all_subsampled_ipsi_radj_cell_gen)

subplot(3,1,3)
histogram(all_contra_radj_gen)

%%  histogram of distributions of model beta values separately
%   for each model term

plt_data = all_betas;
%plt_data = all_radj_gen;

mdl_subset_ind = all_mdl_labels.data_subset == 'contra';
mdl_subset_ind = mdl_subset_ind | all_mdl_labels.data_subset == 'subsampled_ipsi';

region_ind = all_mdl_labels.region == 'ofc';

sig_x1 = all_ps(:, 1) < 0.05;
sig_x2 = all_ps(:, 2) < 0.05;

subset_ind = mdl_subset_ind & region_ind;

%subset_ind = subset_ind & sig_x1 & sig_x2;

figure(5); clf;
num_bins = 20;

[I, C] = findeach( all_mdl_labels, {'region', 'data_subset'}, subset_ind );

axs = plots.panels( numel(I) * size(plt_data, 2) );

for i = 1:numel(I)

cell_ind = I{i};
ncols = size( plt_data, 2 );
term_names = { 'x1', 'x2' };

for j = 1:ncols
    term_name = term_names{j};
    ax = axs((i-1) * ncols + j); hold( ax, 'on' );
    histogram(ax, plt_data(cell_ind,j), num_bins )
    %plot(ax,plt_data(cell_ind,j))
    title( ax, sprintf('%s %s %s', term_name, string(C.region(i)), string(C.data_subset(i))) );
    med = nanmedian( plt_data(cell_ind, j) );
    shared_utils.plot.add_vertical_lines( ax, med );
    text( ax, med, max(get(ax, 'ylim')) - 0.2, sprintf('m = %0.3f', med) );
end

end


shared_utils.plot.match_xlims( axs );

%%  compare median magnitudes of factor (x1 or x2) beta values between contra and ipsi 
% % (or some other combination of datasets)

data_table = all_betas;
%data_table = all_radj_gen;

[I, C] = findeach( all_mdl_labels, {'region'} );
C.p_contra_v_ipsi = nan( numel(I), 1 );
C.med1 = nan( numel(I), 1 );
C.med2 = nan( numel(I), 1 );

contra_subset_ind = all_mdl_labels.data_subset == 'contra';
ipsi_subset_ind = all_mdl_labels.data_subset == 'subsampled_ipsi';

for i = 1:numel(I)
    reg_contra_ind = intersect( find(contra_subset_ind), I{i} );
    reg_ipsi_ind = intersect( find(ipsi_subset_ind), I{i} );

    x1_contra = data_table(reg_contra_ind, 1);
    x1_ipsi = data_table(reg_ipsi_ind, 1);
    C.p_contra_v_ipsi(i) = ranksum( x1_contra, x1_ipsi );
    C.med1(i) = nanmedian( x1_contra );
    C.med2(i) = nanmedian( x1_ipsi );
    
end




%% scatter plot of contra vs ipsi r^2 and contra coeff vs ipsi coeff

figure(2); clf;

term_index = 1;
is_r2 = true;

if ( is_r2 )
    scatter_data = all_radj_gen;
else
    scatter_data = all_betas;
end

[I, C] = findeach( all_mdl_labels, {'region'} );
contra_scatter_ind = all_mdl_labels.data_subset == 'contra';
ipsi_scatter_ind = all_mdl_labels.data_subset == 'subsampled_ipsi';

sig_x = find( all_ps(:, term_index) < 0.05);
%sig_type = 'both';
sig_type = 'any';

term_labels = {'x1', 'x2'};
term_label = term_labels{term_index};

for i = 1:numel(I)
    aj = subplot(2,2,i);
    reg_scatter_contra = intersect( find(contra_scatter_ind), I{i} ); 
    reg_scatter_ipsi = intersect( find(ipsi_scatter_ind), I{i} );

    switch ( sig_type )
        case 'both'
            [~, sig_contra] = intersect( reg_scatter_contra, sig_x );
            [~, sig_ipsi] = intersect( reg_scatter_ipsi, sig_x );
            sig_both = intersect( sig_contra, sig_ipsi );
            reg_scatter_contra = reg_scatter_contra(sig_both);
            reg_scatter_ipsi = reg_scatter_ipsi(sig_both);
        
        case 'either'
            [~, sig_contra] = intersect( reg_scatter_contra, sig_x );
            [~, sig_ipsi] = intersect( reg_scatter_ipsi, sig_x );
            sig_either = union( sig_contra, sig_ipsi );
            reg_scatter_contra = reg_scatter_contra(sig_either);
            reg_scatter_ipsi = reg_scatter_ipsi(sig_either);
        case 'any'

        otherwise
            error( 'Unrecognized sig type');
    end


   contra_scatter{i} = scatter_data(reg_scatter_contra,term_index); %x1 self
   ipsi_scatter{i} = scatter_data(reg_scatter_ipsi,term_index);

   [corrs,pval] = corr(contra_scatter{i},ipsi_scatter{i});

    %contra_scatter{i} = scatter_data(reg_scatter_contra,2); %x2 other
    %ipsi_scatter{i} = scatter_data(reg_scatter_ipsi,2);

    scatter(contra_scatter{i},ipsi_scatter{i},'filled');
    xlabel( 'contra' ); ylabel( 'ipsi' );
    pbaspect([1 1 1])
    if ( is_r2 )
        xlim([-0.25 1.15])
        ylim([-0.2 1.15]);
    else
        xlim([-0.4, 0.4]);
        ylim([-0.4, 0.4]);
    end
    %(R = %0.3f, P = %0.3f)
    title( sprintf( '%s | %s', term_label, char(C.region(i))) );
    %title( sprintf('%s %s (R = %0.3f, P = %0.3f)', term_label, char(C.region(i)), corrs, pval) );
    x=0;
    y=0;
    x=y;
    hold on
   xline(x,'LineWidth',2.5)
   yline(x, 'LineWidth', 2.5)
   plot([0 1],[0 1],'LineWidth',2.5, 'color', 'k')
   plot([-1 0],[-1 0],'LineWidth',2.5, 'color', 'k')
    %labels = {'n(141)','n(364)', 'n(111)', 'n(136)'};
    %lgd = legend(labels);

    total_num_plotted = sum(...
        ~isnan(contra_scatter{i}) & ~isnan(ipsi_scatter{i}));

    S.totalcellsregion_contra(i) = sum(~isnan(contra_scatter{i}))
    S.totalcellsregion_ipsi(i) = sum(~isnan(ipsi_scatter{i}))
    S.totalcellsregion(i) = numel(contra_scatter{i})
    S.totalcellsplotted(i) = total_num_plotted; 

end 

% correlation bw contra and ipsi 
all_regions_corrs = {};
all_region_pval_ipsi_contra_coef = {};

    for i =  1:numel(I)
        [corrs,pval] = corr(contra_scatter{i},ipsi_scatter{i});
        all_regions_corrs{end+1} = corrs;           
        all_region_pval_ipsi_contra_coef{end+1} = pval;
    end 

%%  beta diffs contra - ipsi

do_box = false;
do_hist = false;
term_index = 1;

 contra_ind = all_mdl_labels.data_subset == 'contra';
 ipsi_ind = all_mdl_labels.data_subset == 'subsampled_ipsi';
% 
 diffs = all_betas(contra_ind, :) - all_betas(ipsi_ind, :);
%diffs=mdl_betas(:,1)-mdl_betas(:,2);

abscontra = abs(all_betas(contra_ind, 1));
absipsi = abs(all_betas(ipsi_ind, 1));

differences = abscontra - absipsi;

abs_diffs = abs( diffs ); 

%   stat for median test against 0.
[reg_I, reg_C] = findeach( all_mdl_labels(contra_ind, :), 'region' );
[sr_res, ~, stats] = cellfun( @(x) signrank(differences(x, term_index)), reg_I );
reg_C.sr_ps = sr_res;
%reg_C.sr_z = [stats.zval]';

clf;
if ( do_box )
    [I, id, C] = rowsets( 2, all_mdl_labels(contra_ind, :), 'data_subset', 'region' );
    L = plots.cellstr_join( C );
    [PI, PL] = plots.nest2( id, I, L );
    
    axs = plots.panels( numel(PI) );
    for i = 1:numel(PI)
        plots.box( axs(i), differences(:, term_index), PI{1}, PL{i, 2}, PL{i, 1} );
    end
elseif ( do_hist )
    [I, id, C] = rowsets( 1, all_mdl_labels(contra_ind, :), {'data_subset', 'region'} );
    L = plots.cellstr_join( C );

    axs = plots.panels( numel(I) );
    term_splits = nan( numel(axs), 1 );
    term_split_table = vertcat( C{:} );
    for i = 1:numel(I)
        plots.hist( axs(i), differences(I{i}, term_index), L{i} );
        term_splits(i) = prctile( differences(I{i}, term_index), 80 );
    end
    term_split_table.splits = term_splits;
    shared_utils.plot.match_ylims( axs );
%     shared_utils.plot.match_xlims( axs );

else
    [I, id, C] = rowsets( 3, all_mdl_labels(contra_ind, :), 'data_subset', 'region', {} );
    L = plots.cellstr_join( C );

    axs = plots.simplest_barsets( differences(:, term_index), I, id, L ...
        , 'summary_func', @nanmedian ...
        , 'error_func', @(x) nan ...
    );
    shared_utils.plot.match_ylims( axs );
    title(axs, 'slope differences between contra & ipsi');
     %ylim( [-0.1, 0.1] );
end

%% mean coefficient values

term_index = 1;
bar_data = abs( all_radj_gen(:, term_index) );

mask = all_ps(:, term_index) < 0.05;
mask(:) = true;
    
figure(1); clf;
[I, id, C] = rowsets( 3, all_mdl_labels, {'region'}, {'data_subset'}, {'region'} ...
    , 'mask', mask ...
);
L = plots.cellstr_join( C );
axs = plots.simplest_barsets( bar_data, I, id, L ...
    , 'summary_func', @nanmean ...
    , 'error_func', @(x) plotlabeled.nansem(x) ...
);

%%  plot pie chart showing how many cells have both ipsi and contra x1 or x2

scatter_data = all_betas;

reg_scatter_ind = all_mdl_labels.data_subset == 'contra';
reg_scatter_ind = all_mdl_labels.data_subset == 'subsampled_ipsi';
[I, C] = findeach( all_mdl_labels, {'region', 'data_subset'}, reg_scatter_ind );
% reg_scatter_ind = all_mdl_labels.data_subset == 'subsampled_ipsi';
%subsampled_ipsi_scatter_ind = all_mdl_labels.data_subset == 'subsampled_ipsi';

figure(2); clf;
p_thresh = 0.05;

for i = 1:numel(I)
    aj = subplot(2,2,i);
    reg_scatter_ind = I{i};

    sig_x1 = all_ps(reg_scatter_ind, 1) < p_thresh;
    sig_x2 = all_ps(reg_scatter_ind, 2) < p_thresh;

    both_sig = sum( sig_x1 & sig_x2 );
    x1_sig = sum( sig_x1 & ~sig_x2 );
    x2_sig = sum( sig_x2 & ~sig_x1 );
    rest = sum( ~sig_x1 & ~sig_x2 );

    pie([x1_sig, x2_sig, both_sig, rest]);
    title( sprintf('Dist. of distance monitoring cells (%s | %s)', C.region(i), C.data_subset(i)) );
    labels = {'both x1 & x2','x1 only', 'x2 only', 'neither'};
    lgd = legend(labels);

end 
%%

% to find cells firing across 1 session(day) averaged into different bins

clf;
figure(1);

[I, C] = findeach( spike_labels, 'region' );
psth_by_region = cell( numel(I), 1 );
psth_labels_by_region = cell( numel(I), 1 );
psth_distances = cell( numel(I), 1 );

for i = 1:numel(I)
    
    ax = subplot(4,1,i);
    psth_by_region(i) = {vertcat(allpsth{I{i}})};
    psth_labels_by_region(i) = {vertcat(all_psth_labels{I{i}})};
    psth_distances(i) = {vertcat(all_psth_distances_to_self_origin{I{i}})};

    %close_distance = psth_distances{i} < 100;
    %med_distance = psth_distances{i} >= 10 & psth_distances{i} < 20;

%     close_distance = close_distance & psth_labels_by_regions{i}.id_m1 == 'm1_lynch';
 
    %   compute a sum over elements of allpsth for a given region,
    %   identified by I{i}
    avg_psth_by_region(i) = {mean(psth_by_region{i})};
%     avg_psth_by_region(i) = {mean(psth_by_region{i}(close_distance, :))};

    plot(bin_t, avg_psth_by_region{i});
    title( C.region(i) );

    pos_psth_elements = all_psth_m1_hemifield_delta_pos{i} >0;

    neg_psth_elements = all_psth_m1_hemifield_delta_pos{i} < 0;
 
end


deg_pos_psth = vertcat( all_psth_m1_hemifield_delta_pos{:} );
deg_pos_psth = bfw.px2deg( deg_pos_psth );



%% binning data based on close or distance fixations

clf;
figure(1);

[I, C] = findeach( spike_labels, 'region' );

for i = 1:numel(I)
 ax = subplot(4,1,i);

 close_distance = deg_pos_psth(I{i}) < 5;
 med_distance = deg_pos_psth(I{i}) >= 5 & deg_pos_psth(I{i}) < 10;
 far_distance = deg_pos_psth(I{i}) >=10 & deg_pos_psth(I{i}) < 15;
 furthest_distance = deg_pos_psth(I{i}) >= 15 & deg_pos_psth(I{i}) <20;

 sum_close_distance(i) = sum(deg_pos_psth(I{i}) < 5);
 sum_med_distance(i) = sum(deg_pos_psth(I{i}) >= 5 & deg_pos_psth(I{i}) < 10);
 sum_far_distance(i) = sum(deg_pos_psth(I{i}) >=10 & deg_pos_psth(I{i}) < 15);
 sum_furthest_distance(i) = sum(deg_pos_psth(I{i}) >= 15 & deg_pos_psth(I{i}) <20);

end 


pie([sum_close_distance(1,1), sum_med_distance(1,1), sum_far_distance(1,1), sum_furthest_distance(1,1)])
labels = {'Close Distance','Med Distance', 'Far Distance', 'Furthest Distance'};
lgd = legend(labels);

%pie([sum_close_distance(1,1), sum_med_distance(1,1), sum_far_distance(1,1), sum_furthest_distance(1,i)], '%.2f%%')
%title( C.region(i) );

%% M1's fixations must be 20 degrees from the center of M2's eyes (0 degrees)
% seperate by contra and ipsi



[I, C] = findeach( spike_labels, 'id_m1' );


m1_reject = [];
for i = 1:numel(I)
     ax = subplot(4,1,i);

    m1_fixation_events_regions(i) = {deg_pos_psth(I{i}) < 20};
    sum_m1_fixation_events(i) = sum(m1_fixation_events_regions{i});


   m1_contra_fix(i) = {deg_pos_psth(I{i}) > 0  & deg_pos_psth(I{i}) <20};
   m1_ipsi_fix(i) = {deg_pos_psth(I{i}) < 0  & deg_pos_psth(I{i}) > -20};
   m1_reject(i) = pnz( ~m1_contra_fix{i} & ~m1_ipsi_fix{i} );

   %need to change to be =< 20 and =< -20 so that it includes those that
   %are right after the cusp

   sum_contra_fix(i) = sum(m1_contra_fix{i} > 0);
   sum_ipsi_fix(i) = sum(m1_ipsi_fix{i});


    


   pie([sum_ipsi_fix(i), sum_contra_fix(i)], '%.2f%%');
   labels = {'Ipsilateral','Contralateral'};
   title( ax,  C.id_m1(i) );
   lgd = legend(labels);

end

sum_m1_fixation_events = transpose(sum_m1_fixation_events);
m1_contra_fix = transpose(m1_contra_fix);
m1_ipsi_fix = transpose(m1_ipsi_fix);
sum_contra_fix = transpose(sum_contra_fix);
sum_ipsi_fix = transpose(sum_ipsi_fix);




%%

function [contra_ipsi_label, distance_labels] = label_fixation_distances(hemifield_deg_pos, eye_dist)

m1_contra_fix = hemifield_deg_pos > 0  & hemifield_deg_pos < 20;
m1_ipsi_fix = hemifield_deg_pos < 0  & hemifield_deg_pos > -20;

contra_ipsi_label = cell( size(hemifield_deg_pos) );
contra_ipsi_label(:) = {''};

contra_ipsi_label(m1_contra_fix) = {'contra'};
contra_ipsi_label(m1_ipsi_fix) = {'ipsi'};

distance_windows = { [0, 5], [5, 10], [10,15], [15,20]};
distance_window_labels = { 'close', 'medium', 'far', 'furthest' };

distance_labels = cell( size(contra_ipsi_label) );
distance_labels(:) = {''};

for i = 1:numel(distance_windows)
    dist_win = distance_windows{i};
    within_eye_dist =  eye_dist >= min(dist_win)  & eye_dist <= max(dist_win);
    distance_labels(within_eye_dist) = distance_window_labels(i);
end

end

%% seperating fixations for M2 and applying 20 degree criteria


function [m2_dist, would_pass_sample_crit] = m2_distance_during_m1_fixation(m1_t0, m1_t1, time_file, pos_file, roi_file)

would_pass_sample_crit = false;

if ( isnan(m1_t0) || isnan(m1_t1) )
    m2_dist = nan;
    return
end

m2_pos_subset = get_m2_fixation_position( m1_t0, m1_t1, time_file.t, pos_file.m2 );
m2_pos_centroid = nanmean( m2_pos_subset, 2 );

%   center of m1's eyes from m2's perspective
m1_eye_center = shared_utils.rect.center(roi_file.m2.rects('eyes_nf'));
%   distance (in pixels) for mean m2 gaze position to m1's eyes
m2_dist_to_m1_eyes = norm( m2_pos_centroid(:) - m1_eye_center' );
%   convert to degrees
m2_dist_to_m1_eyes = bfw.px2deg( m2_dist_to_m1_eyes );
m2_pos_subset = bfw.px2deg( m2_pos_subset );

%   check to see if this m2 distance value should be rejected, either
%   because it is outside of 20degrees or because m2 is saccading
m2_is_fixating = samples_within_tolerance_of_center( m2_pos_subset, 1, 0.9 );
% m2_is_within_dist_crit = m2_dist_to_m1_eyes <= 20;
m2_is_within_dist_crit = true;
would_pass_sample_crit = m2_dist_to_m1_eyes <= 20 && m2_is_fixating;

if ( m2_is_fixating && m2_is_within_dist_crit )
    m2_dist = m2_dist_to_m1_eyes;
else
    m2_dist = nan;
end

end

%%  

function contra_ipsi_mask = make_real_contra_ipsi_mask(m1_hemifield_origin_delta_pos, sessions, mdl_prefix)

assert( numel(sessions) == size(m1_hemifield_origin_delta_pos, 1) );

contra_ipsi_mask = false( size(m1_hemifield_origin_delta_pos, 1), 1 );
    
%   If not commented out, then ignore contra vs ipsi and just look at all
%   fixations.
% contra_ipsi_mask(:) = true;

switch ( mdl_prefix )
    case 'contra'
        contra_ipsi_mask = m1_hemifield_origin_delta_pos(:, 1) > 0; 
    case 'ipsi'
       contra_ipsi_mask = m1_hemifield_origin_delta_pos(:, 1) < 0;
       case 'all_contra_ipsi'
        contra_ipsi_mask = m1_hemifield_origin_delta_pos(:, 1);
    case 'subsampled_ipsi'
    otherwise
        error( 'Unrecognized mdl prefix: "%s".', mdl_prefix );
end    

if ( strcmp(mdl_prefix, 'subsampled_ipsi') )
    I = findeach( sessions, 1 );
    for i = 1:numel(I)
        is_ipsi = m1_hemifield_origin_delta_pos(I{i}, 1) < 0;
        num_contra = sum( m1_hemifield_origin_delta_pos(I{i}, 1) > 0 );
        num_ipsi = sum( is_ipsi );
        ipsi_indices = I{i}(is_ipsi);
        if num_contra > num_ipsi
            keep_ipsi_indicies = ipsi_indices;
        else 
            keep_ipsi_indicies = randsample(ipsi_indices,num_contra);
        end 
        contra_ipsi_mask(keep_ipsi_indicies) = true;
    
        %   hint: `help randsample`
        %   define the `kept_ipsi_indices` by sub-selecting from `ipsi_indices`
        %   contra_ipsi_mask(kept_ipsi_indices) = true;
    end
end


end
