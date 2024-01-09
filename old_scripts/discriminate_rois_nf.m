%data_p = '~/Downloads';

data_p = fullfile(eisg.util.project_path, 'processed_data');

sorted = shared_utils.io.fload( fullfile(data_p, 'sorted_neural_data_social_gaze.mat') );
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );
ct_labels = shared_utils.io.fload( fullfile(data_p, 'celltype_labels_p2v_combined.mat') );

%%

[spike_ts, ~, spike_labels] = linearize_sorted( sorted );
apply_fcat_celltype_labels( spike_labels, ct_labels );
% spike_labels = apply_cell_type_labels(...
%   spike_labels, ct_labels.label_mat, ct_labels.label_mat_cols );

events = add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

%%  psth

rois = { 'whole_face', 'right_nonsocial_object_whole_face_matched' };

evt_mask = find( events.labels, [{'m1'}, rois] );
spk_mask = find( spike_labels, {'valid-unit', 'maybe-valid-unit'} );

min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;

[psth_matrix, psth_labels, t] = compute_psth(...
    spike_ts, spike_labels, spk_mask ...
  , evts, events.labels, evt_mask ...
  , min_t, max_t, bin_width );



%%  discriminate between pairs of rois

over_time = true;

if ( over_time )
  if ( 0 )
    bin_step = 0.1;
    bin_size = 0.15;
    t_wins = -0.5:bin_step:0.5-bin_size;
    t_wins = arrayfun( @(x) [x, x+bin_size], t_wins, 'un', 0 );
  else
    t_wins = arrayfun( @(x, y) [x, y], t(1:end-1), t(2:end), 'un', 0 );
  end
else
  t_wins = { [0, 0.5] };
end

fit_res_over_time = cell( size(t_wins) );
ps_over_time = cell( size(t_wins) );

for i = 1:numel(t_wins)
 
t_win = t_wins{i};

t_mask = t >= t_win(1) & t < t_win(2);
t_mean = mean( psth_matrix(:, t_mask), 2 );
t_subset = psth_matrix(:, t_mask);
% t_subset = psth_matrix;

discrim_x = t_mean; % use single value for each observation
% discrim_x = t_subset; % use vector of values for each observation

rois = combs( psth_labels, 'roi' );
roi_pairs = bfw.pair_combination_indices( numel(rois) );

as = arrayfun( @(x) rois{x}, roi_pairs(:, 1), 'un', 0 );
bs = arrayfun( @(x) rois{x}, roi_pairs(:, 2), 'un', 0 );

% discriminate between a and b, for each a in as and b0 in bs
each_I = findall( psth_labels, {'uuid', 'looks_by'} );
[fit_res, fit_labels] = binary_fitcdiscrs( ...
  discrim_x, psth_labels, each_I, as, bs, 'roi' );

fit_res_over_time{i} = fit_res;
ps_over_time{i} = cellfun( @(x) x.p, fit_res );

end
%%
save( fullfile(data_p, 'roi_classification_over_time.mat'), 'fit_res_over_time', 'ps_over_time', 'fit_labels');

%%  plot single unit psth, overlaying significantly decoded bins

ps = horzcat( ps_over_time{:} );
is_sig = double( ps < 0.05 );
assert_ispair( ps, fit_labels );

[fig_I, fig_C] = findall( psth_labels, {'uuid', 'region'} );
% fig_I = fig_I(1);

for i = 1:numel(fig_I)
  fi = fig_I{i};
  
  pcats = {'cell-type', 'looks_by', 'region', 'uuid' };
  [p_I, p_C] = findall( psth_labels, pcats, fi );
  
  shp = plotlabeled.get_subplot_shape( numel(p_I) );
  for j = 1:numel(p_I)
    ax = subplot( shp(1), shp(2), j );
    cla( ax ); hold( ax, 'on' );
    
    gcats = { 'roi' };
    [g_I, g_C] = findall( psth_labels, gcats, p_I{j} );
    psth_means = bfw.row_nanmean( psth_matrix, g_I );
    
    hs = gobjects( numel(g_I), 1 );
    for k = 1:size(psth_means, 1)
      hs(k) = plot( gca, t, psth_means(k, :) );
    end
    
    match_sig = find( fit_labels, fig_C(:, i) );
    if ( numel(match_sig) == 1 )
      sig_sub = find( is_sig(match_sig, :) );
      plot( gca, t(sig_sub), repmat(max(get(gca, 'ylim')), numel(sig_sub), 1), 'k*' );
    end
    
    rep_join = @(s) strrep( fcat.strjoin(s, ' | '), '_', ' ' );
    
    legend( hs, rep_join(g_C) );
    title( rep_join(p_C(:, j)) );
  end
  
  if ( 1 )
    reg_subdir = fig_C{2, i};
    save_labs = prune( psth_labels(fi) );
    save_p = fullfile( eisg.util.project_path(), 'data/plots/decoding_over_time_psth', dsp3.datedir );
    save_p = fullfile( save_p, reg_subdir );
    dsp3.req_savefig( gcf, save_p, save_labs, 'uuid' );
  end
end

%% plot probability time course

ps = horzcat( ps_over_time{:} );
is_sig = double( ps < 0.05 );
assert_ispair( ps, fit_labels );

% proportions of significant cells per region and roi
mask = findnone( fit_labels, {'<cell-type>', 'b'} );
[prop_labs, prop_I] = retaineach(...
  fit_labels, {'region', 'looks_by', 'roi', 'cell-type'}, mask );
props = matrix_proportions( is_sig, prop_I );

pl = plotlabeled.make_common();
pl.x = t(1:end-1);
axs = pl.lines( props, prop_labs, 'cell-type', {'roi', 'looks_by', 'region'} );
ylabel( axs(1), 'Prop Significant a vs b.' );


%%  plot results

% extract p value from classification
ps = cellfun( @(x) x.p, fit_res );
is_sig = ps < 0.05;

% proportions of significant cells per region and roi
mask = findnone( fit_labels, {'<cell-type>'} );
[prop_labs, prop_I] = retaineach(...
  fit_labels, {'region', 'looks_by', 'roi', 'cell-type'}, mask );
props = proportions( is_sig, prop_I );

pl = plotlabeled.make_common();
axs = pl.bar( props, prop_labs, {}, 'cell-type', {'roi', 'looks_by', 'region'} );
ylabel( axs(1), 'Prop Significant a vs b.' );