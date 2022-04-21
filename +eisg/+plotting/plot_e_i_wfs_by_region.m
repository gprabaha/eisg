function plot_e_i_wfs_by_region(sorted_neural_data, varargin)

defaults = eisg.util.make_analysis_params_struct();

params = shared_utils.general.parsestruct( defaults, varargin );

regions = unique( {sorted_neural_data.region} );

growing_wfs = [];
growing_labels = [];

for region = regions
    params.region = char(region);
    ei_labels = eisg.util.fetch_ei_labels( region, params );
    label_mat = ei_labels.label_mat;
    label_cols = ei_labels.label_mat_cols;
    uuid_col_in_ei_mat = strcmp( label_cols, 'uuid' );
    file_ind_col_in_ei_mat = strcmp( label_cols, 'file_ind' );
    unit_ind_col_in_ei_mat = strcmp( label_cols, 'unit_ind' );
    ei_label_col_in_ei_mat = strcmp( label_cols, 'ei_label' );
    regional_uuids = label_mat(:, uuid_col_in_ei_mat)';
    for uuid_index = 1:length(regional_uuids)
        uuid_file_ind = double( string( label_mat(uuid_index, file_ind_col_in_ei_mat) ) );
        uuid_unit_ind = double( string( label_mat(uuid_index, unit_ind_col_in_ei_mat) ) );
        uuid_wf = sorted_neural_data(uuid_file_ind).templates(uuid_unit_ind,:);
        uuid_timepoints = sorted_neural_data(uuid_file_ind).template_timepoints;
        uuid_ei_label = label_mat(uuid_index, ei_label_col_in_ei_mat);
        uuid_labels = [region, label_mat(uuid_index, uuid_col_in_ei_mat), uuid_ei_label];
        growing_wfs = [growing_wfs; uuid_wf];
        growing_labels = [growing_labels; uuid_labels];
    end
end
wf_label_categories = {'region', 'uuid', 'ei_label'};
wf_labels = fcat.from( growing_labels, wf_label_categories );

normalized_wfs = normalize_wfs( growing_wfs, wf_labels );

% Plotting

% For the average traces
pl = plotlabeled.make_common();
x_vec = linspace(-3000, 3000, 241);
pl.x = x_vec(2:end);
pl.summary_func = @nanmean;
pl.error_func = @plotlabeled.nansem;
pl.main_line_width = 2;
pl.error_line_width = 1;
[axs, hs, ids] = pl.lines( normalized_wfs, wf_labels, 'ei_label', 'region' );

% For individual wfs
wf_alpha = 0.15;
set(findobj(gcf, 'type', 'legend'), 'autoupdate', false);
for i = 1:numel(ids)
  % for each panel
  ax = axs(i);
  hold( ax, 'on' );
  indices = ids{i};
  for j = 1:numel(indices)
      % for each line
      plot_line = hs{i}(j);
      trace_color =  get(plot_line, 'color');
      trace_color = [trace_color wf_alpha];
      index = indices{j};
      subset = normalized_wfs(index, :);
      h = plot( ax, x_vec(2:end), subset' );
      set( h, 'linewidth', 0.5 );
      set( h, 'color', trace_color );
  end
end

save_p = '~/Desktop/example_plots';
dsp3.req_savefig( gcf, save_p, wf_labels, {}, 'social_gaze_e_i_efs');

mn = min( normalized_wfs(:) );
mx = max( normalized_wfs(:) );
shared_utils.plot.set_ylims( axs, [mn, mx] );

% pcats = { 'region' };
% gcats = { 'ei_label' };
% 
% x_vec = 1:size( normalized_wfs, 2 );
% 
% panel_I = findall( wf_labels, pcats );
% shp = plotlabeled.get_subplot_shape( numel(panel_I) );
% 
% for i = 1:numel(panel_I)
%     % for each panel
%     ax = subplot( shp(1), shp(2), i );
%     cla( ax );
%     hold( ax, 'on' );
% 
%     % plot a line for each group.
%     [group_I, group_C] = findall( wf_labels, gcats, panel_I{i} );
%     colors = hsv( numel(group_I) );
% 
%     for j = 1:numel(group_I)
%         % for each group
%         gi = group_I{j};
%         group_label = strjoin( group_C(:, j), ' | ' );
% 
%         wfs = normalized_wfs(gi, :);
%         h_wfs = plot( ax, x_vec, wfs' );
% 
%         mean_wfs = nanmean( wfs, 1 );
%         err_wfs = nanstd( wfs, [], 1 );
% 
%         hs = gobjects( 0 );
%         hs(end+1) = plot( ax, x_vec, mean_wfs, 'displayname', group_label );
%         hs(end+1) = plot( ax, x_vec, mean_wfs + err_wfs );
%         hs(end+1) = plot( ax, x_vec, mean_wfs + err_wfs );
% 
%         set( hs, 'color', colors(j, :) );
%         set( h_wfs, 'color', colors(j, :) );    
%     end
% end

end

function normalized_wfs = normalize_wfs(growing_wfs, wf_labels)

normalized_wfs = growing_wfs;

regions = wf_labels('region');
for region = regions'
    regional_wf_inds = find(wf_labels, region);
    regional_wfs = growing_wfs(regional_wf_inds, :);
    regional_wfs_norm = regional_wfs./max( abs(regional_wfs), [], 2 );
    normalized_wfs(regional_wf_inds, :) = regional_wfs_norm;
end

end