%%

data_dir = '/Users/shadeeleazer/data/brains/intermediates';

pos_files = shared_utils.io.findmat( fullfile(data_dir, 'aligned_raw_samples/position') );

matched_files = bfw.matched_files( pos_files ...
    , fullfile(data_dir, 'aligned_raw_samples/time') ...
    , fullfile(data_dir, 'events') ...
    , fullfile(data_dir, 'rois') ...
    , fullfile(data_dir, 'single_origin_offsets') ... 
    , fullfile(data_dir, 'meta'));

%%
% create loop to create matrix for distances for each run

all_distances = {};
all_delta_pos = {};
all_labels = {};
all_fixation_ts = {};

%   @TODO: Make another aggregate array for the fixation time stamps and
%   append to it, like we did for the distances

%%

%includes elements for each run
for i = 1:length(matched_files)
    %%
    fprintf( '\nProcessing: (%d of %d)', i, length(matched_files));

    one_file = matched_files(i,:);

    pos_file = shared_utils.io.fload( one_file{1} );
    time_file = shared_utils.io.fload( one_file{2} );
    events_file = shared_utils.io.fload( one_file{3} );
    roi_file = shared_utils.io.fload( one_file{4} );
    offsets_file = shared_utils.io.fload( one_file{5} );
    meta_file = shared_utils.io.fload( one_file{6} );

    if ( isfield(pos_file, 'm2')) 
        [dists, eye_dists, delta_pos, roi_labels] = ...
            get_hemifield_center_distance( pos_file, events_file, roi_file, offsets_file, false, 'm1' );

        fix_ts = get_fixation_times( events_file );
    
        %   append to arrays
        all_distances{end+1} = dists;
        all_delta_pos{end+1} = delta_pos;
        all_labels{end+1} = repmat( struct2cell( meta_file )', numel(dists), 1 );
        all_fixation_ts(end+1) = { fix_ts };
    end
end


%%

deg_pos = cell( size(all_delta_pos) );
for i = 1:numel(all_delta_pos)
    deg_pos{i} = bfw.px2deg( all_delta_pos{i} );
end

mean_x = cellfun( @(x) nanmean(x(:, 1)), deg_pos );
above_thresh = mean_x > 30;

%%  create matrix of labels by vertically concenating sub-matrices

mat_labs = vertcat( all_labels{:} );
[id_m1, id_m2] = extract_m1_m2_ids( mat_labs(:, 3) );

%%  convert pixel distances to degrees

% deg_pos = cell( size(all_delta_pos) );
% for i = 1:numel(all_delta_pos)
%     deg_pos{i} = bfw.px2deg( all_delta_pos{i} );
% end

deg_pos = vertcat( all_delta_pos{:} );
deg_pos = bfw.px2deg( deg_pos );

%%  plot histogram of left vs right fixations

clf;

nbins = round(sqrt(numel(deg_pos(:,1))));

xele_deg_pos = deg_pos(:, 1);

histogram(xele_deg_pos,nbins);

%   add vertical line for median x degree offset from center
hold on;
shared_utils.plot.add_vertical_lines( gca, nanmedian(deg_pos(:, 1)) );

xlim( [-10, 10] );

figure(1);

%%  split histogram by m1 id

%   example
x = rand( 8, 1 );   %   values
ind = [3, 4, 5];    %   indices into `x`
y = x(ind);         %   subset of `values` selected by `ind`

%%  conditionally flip left vs right so that + means ipsilateral and - means contralateral

[I, C] = findeach( id_m1, 1 );

contra_ipsi_deg_pos = xele_deg_pos;
for i = 1:numel(I)
   if ( strcmp(C{i}, 'm1_lynch') )
        %   flip left and right for lynch, because his chamber's on the
        %   left. add a negaitve sign to flip across and make median
        %   positive
      contra_ipsi_deg_pos(I{i}) = contra_ipsi_deg_pos(I{i});

    elseif ( strcmp(C{i}, 'm1_kuro') )
        %Kuro's chamber is also on the right
        %kuro's sign needed to be flipped bc their chambers are on opposite
        %sides. So kuro's contra is left while lynch's is right
       contra_ipsi_deg_pos(I{i}) = -contra_ipsi_deg_pos(I{i});
           ... 

  else
       error( 'Expected lynch or kuro.' );
   end
end


%can flip between ipsi and contra (when using ipsi is on the right)
%%

%figure(2);
clf;
[I, C] = findeach( id_m1, 1 );

plot_src = contra_ipsi_deg_pos; %   contra vs ipsi
% plot_src = xele_deg_pos;    %   left vs right

for i = 1:numel(I)
    ax = subplot( 2, 1, i );

    edges = -10:1/2:10;
%     edges = [-1, 0, 1];
    
%     histogram(xele_deg_pos(I{i}),round(sqrt(numel(I{i}))));
    histogram(plot_src(I{i}),edges);
    [~, ~, bin_indices] = histcounts(plot_src(I{i}),edges);

    hold on;

    ylim( [0, 500] );

    med = nanmedian(plot_src(I{i}));
    shared_utils.plot.add_vertical_lines( gca, med );

    is_significant = signrank(plot_src(I{i})) < 0.05

%is_significant = true;
    %   determine whether median difference from 0 is significant
    %

    plot_str = sprintf('Median = %0.3f', med);
    if ( is_significant )
        plot_str = sprintf( '%s (*)', plot_str );
    end

    text( ax, med, 490, plot_str );
    
    title( ax, strrep(C{i}, '_', ' ') );
    xlabel('degrees')
    ylabel('frequency')

%         histogram(deg_pos(I{i}, 1),round(sqrt(numel(I{i}))));
     

%     bx = subplot( 2, 1, i );
%       histogram(xele_deg_pos(I{i}),round(sqrt(numel(I{i}))));


    %   use the index I{i} to select elements of `deg_pos` to plot a
    %   histogram over



end

%figure(1);

%%  compute medians and stat per day

table_labs = array2table( ...
    [mat_labs, id_m1], 'VariableNames' ...
    , {'unified_filename', 'date', 'session', 'run', 'task_type', 'run_number', 'id_m1'} );

% [I, C] = findeach( table_labs, 'session' );
[day_labels, I] = retaineach( table_labs, 'session' );

day_medians = nan( numel(I), 1 );
median_is_sig = false( numel(I), 1 );

for i = 1:numel(I)
  day_medians(i) = nanmedian(plot_src(I{i})); 
  % use signrank to determine whether distribution has median different
  % from 0.
  median_is_sig(i) = signrank(plot_src(I{i})) < 0.05;
  median_pvalues(i) = signrank(plot_src(I{i}))';
end

%%  plot percentages of positively vs negatively biased days, per animal (m1)

[monk_labels, I] = retaineach( day_labels, 'id_m1' );

figure(1);
for i = 1:numel(I)
    ax = subplot( 2, 1, i );
    %   using the `I` variable to select rows for a given animal,
    %   determine whether the median distances are positive or negative,
    %   then count the number of positives and negatives and plot a pie
    %   chart
    pos_neg_med = (day_medians(I{i}) < 0);
    neg_medians = sum(day_medians(I{i}) < 0);
    positive_medians = sum(day_medians(I{i}) > 0);

   pie([neg_medians, positive_medians], '%.2f%%');
   labels = {'Ipsilateral','Contralateral'};
   title( ax, strrep(monk_labels.id_m1{i},'_', ' ') );
   lgd = legend(labels);
 
end

% for Lynch on most days there is a bias towards the ipsi side. 4 neagtive
% means (bias towards contra side). 10 days there is a ipsi bias
%for Kuro 52 days there is an ipsi bias and 1 day there is a contra bias
%%  plot percentages of significantly vs not significantly biased days, per animal (m1)

[monk_labels, I] = retaineach( day_labels, 'id_m1' );

figure(2);
for i = 1:numel(I)
    ax = subplot( 2, 1, i );
    %   using the `I` variable to select rows for a given animal,
    %   count the number of significant medians and plot a pie
    %   chart.
    %   sum of median_is_sig will give you the number of significant
    %   medians.
    sig_positive_med = sum(median_is_sig(I{i}) > 0);
    sig_neg_med = sum(median_is_sig(I{i}) == 0);


   pie([sig_neg_med, sig_positive_med], '%.2f%%');
   labels = {'Insignificant bias', 'Significant bias'};
   title( ax, strrep(monk_labels.id_m1{i}, '_','' ));
   lgd = legend(labels);

end
%% Using p-values generated from signrank test to determine significant, 
% insignificant and not significant
[monk_labels, I] = retaineach( day_labels, 'id_m1' );

figure(3);

for i = 1:numel(I)
    ax = subplot( 2, 1, i );
    
    pvaluesign_notsig = (median_pvalues(I{i}) < 0.05)';
    pvaluesign = sum(median_pvalues(I{i}) < 0.05);
    pvaluenotsign = sum(median_pvalues(I{i}) > 0.05);

    pie([pvaluesign, pvaluenotsign]);
    labels = {'Significant bias', 'Not Significant bias'};
    title( ax, strrep(monk_labels.id_m1{i}, '_','' ));
    lgd = legend(labels);
end


%% Seperate significant days into ipsi or contra lateral
[monk_labels, I] = retaineach( day_labels, 'id_m1' );

figure(4);

for i = 1:numel(I)
    ax = subplot (2,1,i);
    
    sig_ipsi = sum((day_medians(I{i}) < 0) & (median_pvalues(I{i}) < 0.05)');
    sig_contra = sum((day_medians(I{i}) > 0) & (median_pvalues(I{i}) < 0.05)');

  
    pie([sig_ipsi, sig_contra, pvaluenotsign]);
    labels = {'Significant Ipsilateral Bias', 'Significant Contralateral Bias', 'Not Significant Bias'};
    title(ax, strrep(monk_labels.id_m1{i}, '_','' ));
    lgd = legend(labels);
end

  


%%
all_distances = shared_utils.io.fload( 'All_runs_m1_eye_fixation_distances.mat' );

%   for each run, compute a mean fixation distance for that run.
%   over runs, compute a mean of the means
%   visualize this mean with a bar plot

all_means= [];

for i = 1:length(all_distances)
    means = mean(all_distances{i}, 'omitnan');
    all_means(end+1) = means;
end

mean_of_means = mean(all_means);

%bar( mean_of_means );

bar(mean_of_means)

%% Find the distribution and determine if there is any biased in the fixation of each animal



