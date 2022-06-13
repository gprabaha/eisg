%%

data_p = '/Users/prabaha/repositories/eisg/processed_data';

sorted = shared_utils.io.fload( fullfile(data_p, 'sorted_neural_data_social_gaze.mat') );
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );

%%
[unit_spike_ts, unit_wfs, spike_labels] = linearize_sorted( sorted );
bfw.add_monk_labels( spike_labels );

%%
[unit_wf_features, feature_labels] = get_wf_features_from_sorted( sorted );
bfw.add_monk_labels( feature_labels );

%%
events = add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

%%

% get the mean firing activity of all units here. then use it as a feature
% for the k-means classification

rois_all_fix = { 'everywhere', 'right_nonsocial_object_whole_face_matched', 'whole_face' };
evt_mask = find( events.labels, [{'m1'}, rois_all_fix] );

spk_mask = rowmask( spike_labels );

min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;

[mean_fr_all_fix, mean_fr_labs_all_fix] = baseline_psth_stats( ...
  evts, events.labels, evt_mask, unit_spike_ts, spike_labels, spk_mask ...
  , 'min_t', min_t ...
  , 'max_t', max_t ...
  , 'bin_width', bin_width ...
);

non_roi_fix = { 'everywhere' };
evt_mask = find( events.labels, [{'m1'}, non_roi_fix] );

[mean_fr_non_roi_fix, mean_fr_labs_non_roi_fix] = baseline_psth_stats( ...
  evts, events.labels, evt_mask, unit_spike_ts, spike_labels, spk_mask ...
  , 'min_t', min_t ...
  , 'max_t', max_t ...
  , 'bin_width', bin_width ...
);

%%

% [unit_wf_features_with_all_fix_fr, feature_labels_with_all_fix_fr] = ...
%   append_mean_fr_in_wf_features( unit_wf_features, feature_labels, mean_fr_all_fix, mean_fr_labs_all_fix );
% 
% [unit_wf_features_with_non_roi_fix_fr, feature_labels_with_non_roi_fix_fr] = ...
%   append_mean_fr_in_wf_features( unit_wf_features, feature_labels, mean_fr_non_roi_fix, mean_fr_labs_non_roi_fix );

%% Classify celltypes based on peak-to-valley, by regions

n_cluster_range = [2 5];

% 'b' stands for 'broad, 'm' for 'medium', and 'n' for 'narrow' spiking
% units

% specify a list of features here
% feature_list = {'peak_to_valley', 'halfwidth'};
list_of_all_features = {
    'valley_ind', ...
    'peak_ind', ...
    'peak_to_valley', ...
    'halfwidth', ...
    'peak_trough_ratio', ...
    'recovery_slope', ...
    'repolarization_slope'
};
feature_list = {'peak_to_valley'};
validity_filter = {'valid-unit', 'maybe-valid-unit'};

% get the specified features for each region
num_features = numel( feature_list );
features_mask = find( feature_labels, [feature_list, validity_filter] );
[regional_feature_inds, regions] = findall( feature_labels, 'region', features_mask );
% subset_inds = regional_feature_inds;
subset_inds = { regional_feature_inds{1}; vertcat(regional_feature_inds{2:4}) };

% get indices for p2v values for the wfs
p2v_inds = cellfun( @(x) find(feature_labels, 'peak_to_valley', x), subset_inds, 'un', 0 );

% use pca or not | default is false
use_pca = true;

% do the k-means celltype classification
[celltype_labels, feature_mat]= eisg.celltype_class.do_feature_based_wf_class_for_each_subset(...
  unit_wf_features, feature_labels, subset_inds, num_features, p2v_inds, n_cluster_range, use_pca);

save(fullfile( [ data_p 'celltype_labels_p2v_combined.mat']), celltype_labels);

%% Plot wfs of celltypes across regions

% fig saving parameters
save_p = fullfile( eisg.util.project_path, 'data/plots/celltype_wfs', dsp3.datedir );
do_save = false;

% add cell-type labels to a temp spike_labels variable
[uuid_I, uuid_C] = findall( celltype_labels, {'uuid', 'cell-type'} );
temp_fr_labels = spike_labels;
match_inds = bfw.find_combinations( temp_fr_labels, uuid_C(1, :) );
addcat( temp_fr_labels, 'cell-type' );
for i = 1:numel(match_inds)
  setcat( temp_fr_labels, 'cell-type', uuid_C{2, i}, match_inds{i} );
end

mask = pipe( rowmask(spike_labels) ...
  , @(m) find(spike_labels, validity_filter, m) ...
  , @(m) findnone(spike_labels, '<cell-type>', m) ...
);

gcats = {'cell-type' };
% pcats = {'region', 'id_m1'};
pcats = {'region'};

fig_I = findall_or_one( spike_labels, {}, mask);
for i = 1:numel(fig_I)
  fi = fig_I{i};
  
  [panel_I, panel_C] = findall( spike_labels, pcats, fi );
  axs = gobjects( numel(panel_I), 1 );
  shape = plotlabeled.get_subplot_shape( numel(panel_I) );
  
  match_group_C = combs( spike_labels, gcats, fi );
  match_group_C = fcat.strjoin( match_group_C );
  colors = hsv( numel(match_group_C) );
  colors(:,4) = 0.3;
  
  for j = 1:numel(panel_I)
    axs(j) = subplot( shape(1), shape(2), j );
    
    [group_I, group_C] = findall( spike_labels, gcats, panel_I{j} );
    query_group_C = fcat.strjoin( group_C );
    [~, ind] = ismember( query_group_C, match_group_C );
    
    elements = cellfun( @(x) unit_wfs(x, :), group_I, 'un', 0 );

    hs = gobjects( numel(elements), 1 );
    for k = 1:numel(elements)
      h = plot( axs(j), 1:size(elements{k}, 2), elements{k} );
      hold( axs(j), 'on' );
      set( h, 'color', colors(ind(k), :) );
      hs(k) = h(1);
    end
    
    legend( hs, fcat.strjoin(group_C, ' | ') );
    title( axs(j), strjoin(panel_C(:, j), ' | ') );
  end
end
suptitle(['Classification by: ' get_feature_string_from_labels(feature_list)]);
if ( do_save )
  dsp3.req_savefig( gcf, save_p, prune(spike_labels(fi)), 'region', ['celltypes_classified_by' is_pca_string( use_pca ) get_feature_string_from_labels(feature_list)] );
end

%% Plot violinplots comparing mean fr (all fix) of celltypes across regions

% fig saving parameters
save_p = fullfile( eisg.util.project_path, 'data/plots/celltype_mean_fr', dsp3.datedir );
do_save = true;

% add cell-type labels to a temp spike_labels variable
[uuid_I, uuid_C] = findall( celltype_labels, {'uuid', 'cell-type'} );
temp_fr_labels = mean_fr_labs_all_fix;
match_inds = bfw.find_combinations( temp_fr_labels, uuid_C(1, :) );
addcat( temp_fr_labels, 'cell-type' );
for i = 1:numel(match_inds)
  setcat( temp_fr_labels, 'cell-type', uuid_C{2, i}, match_inds{i} );
end

mask = pipe( rowmask(temp_fr_labels) ...
  , @(m) find(temp_fr_labels, validity_filter, m) ...
  , @(m) findnone(temp_fr_labels, '<cell-type>', m) ...
);

pl = plotlabeled.make_common();
axs = pl.violinplot( mean_fr_all_fix(mask, 1), prune( temp_fr_labels(mask) ), 'cell-type', 'region' );
if ( do_save )
  dsp3.req_savefig( gcf, save_p, prune(temp_fr_labels), 'region', ['mean_fr_all_fix_across_regional_celltypes_classified_by' is_pca_string( use_pca ) get_feature_string_from_labels(feature_list)] );
end

%% Plot violinplots comparing mean fr (non roi fix) of celltypes across regions

% fig saving parameters
save_p = fullfile( eisg.util.project_path, 'data/plots/celltype_mean_fr', dsp3.datedir );
do_save = true;

% add cell-type labels to a temp spike_labels variable
[uuid_I, uuid_C] = findall( celltype_labels, {'uuid', 'cell-type'} );
temp_fr_labels = mean_fr_labs_non_roi_fix;
match_inds = bfw.find_combinations( temp_fr_labels, uuid_C(1, :) );
addcat( temp_fr_labels, 'cell-type' );
for i = 1:numel(match_inds)
  setcat( temp_fr_labels, 'cell-type', uuid_C{2, i}, match_inds{i} );
end

mask = pipe( rowmask(temp_fr_labels) ...
  , @(m) find(temp_fr_labels, validity_filter, m) ...
  , @(m) findnone(temp_fr_labels, '<cell-type>', m) ...
  , @(m) findnone(temp_fr_labels, '', m) ...
);

pl = plotlabeled.make_common();
axs = pl.violinplot( mean_fr_non_roi_fix(mask, 1), prune(temp_fr_labels(mask)), {'cell-type'}, 'region' );
if ( do_save )
  dsp3.req_savefig( gcf, save_p, prune(temp_fr_labels), 'region', ['mean_fr_non_roi_fix_across_regional_celltypes_classified_by' is_pca_string( use_pca ) get_feature_string_from_labels(feature_list)] );
end

%% Classify celltypes based on peak-to-valley, by combined across PFC regions

n_cluster_range = [2 5];

% 'b' stands for 'broad, 'm' for 'medium', and 'n' for 'narrow' spiking
% units

% specify a list of features here
feature_list = {'peak_to_valley', 'halfwidth'};
validity_filter = {'valid-unit', 'maybe-valid-unit'};

% get the specified features for each region
% following lines have been executed before and commented

% num_features = numel( feature_list );
% features_mask = find( feature_labels, [feature_list, validity_filter] );
% [regional_feature_inds, regions] = findall( feature_labels, 'region', features_mask );

% collapse all PFC inds together
subset_inds = { regional_feature_inds{1}; vertcat(regional_feature_inds{2:4}) };

% get indices for p2v values for the wfs
p2v_inds = cellfun( @(x) find(feature_labels, 'peak_to_valley', x), subset_inds, 'un', 0 );

% use pca or not | default is false
use_pca = false;

% do the k-means celltype classification
[celltype_labels_2, feature_mat]= eisg.celltype_class.do_feature_based_wf_class_for_each_subset(...
  unit_wf_features, feature_labels, subset_inds, num_features, p2v_inds, n_cluster_range, use_pca);

% part of script so collapse regions together; has to be shifted below


%%

% fig saving parameters
save_p = fullfile( eisg.util.project_path, 'data/plots/celltype_wfs', dsp3.datedir );
do_save = true;

% add cell-type labels to a temp spike_labels variable
[uuid_I, uuid_C] = findall( celltype_labels_2, {'uuid', 'cell-type'} );
temp_fr_labels = spike_labels;
match_inds = bfw.find_combinations( temp_fr_labels, uuid_C(1, :) );
addcat( temp_fr_labels, 'cell-type' );
for i = 1:numel(match_inds)
  setcat( temp_fr_labels, 'cell-type', uuid_C{2, i}, match_inds{i} );
end

mask = pipe( rowmask(spike_labels) ...
  , @(m) find(spike_labels, validity_filter, m) ...
  , @(m) findnone(spike_labels, '<cell-type>', m) ...
);

gcats = {'cell-type' };
% pcats = {'region', 'id_m1'};
pcats = {'region'};

fig_I = findall_or_one( spike_labels, {}, mask);
for i = 1:numel(fig_I)
  fi = fig_I{i};
  
  [panel_I, panel_C] = findall( spike_labels, pcats, fi );
  axs = gobjects( numel(panel_I), 1 );
  shape = plotlabeled.get_subplot_shape( numel(panel_I) );
  
  match_group_C = combs( spike_labels, gcats, fi );
  match_group_C = fcat.strjoin( match_group_C );
  colors = hsv( numel(match_group_C) );
  colors(:,4) = 0.3;
  
  for j = 1:numel(panel_I)
    axs(j) = subplot( shape(1), shape(2), j );
    
    [group_I, group_C] = findall( spike_labels, gcats, panel_I{j} );
    query_group_C = fcat.strjoin( group_C );
    [~, ind] = ismember( query_group_C, match_group_C );
    
    elements = cellfun( @(x) unit_wfs(x, :), group_I, 'un', 0 );

    hs = gobjects( numel(elements), 1 );
    for k = 1:numel(elements)
      h = plot( axs(j), 1:size(elements{k}, 2), elements{k} );
      hold( axs(j), 'on' );
      set( h, 'color', colors(ind(k), :) );
      hs(k) = h(1);
    end
    
    legend( hs, fcat.strjoin(group_C, ' | ') );
    title( axs(j), strjoin(panel_C(:, j), ' | ') );
  end
end
suptitle(['Classification by: ' get_feature_string_from_labels(feature_list)]);
if ( do_save )
  dsp3.req_savefig( gcf, save_p, prune(spike_labels(fi)), 'region', ['celltypes_combined_classification_by' is_pca_string( use_pca ) get_feature_string_from_labels(feature_list)] );
end

%%

% region = 'bla';
% subplot(1,2,1);
% valid_mask = find( spike_labels, {region, 'valid-unit', 'maybe-valid-unit'} );
% plot( unit_wfs(valid_mask, :)' );
% title('combined units');
% subplot(1,2,2);
% maybe_valid_mask = find( spike_labels, {region, 'maybe-valid-unit'} );
% plot( unit_wfs(maybe_valid_mask, :)' );
% title('maybe valid units');

%%

unit_mask = find( spike_labels, validity_filter );
[regional_unit_inds, regions] = findall( spike_labels, 'region', unit_mask );
subset_inds = { regional_unit_inds{1}; vertcat( regional_unit_inds{2:4} ) };

features_mask = find( feature_labels, validity_filter );
[regional_feature_inds, regions] = findall( feature_labels, 'region', features_mask );
feature_subset_inds = { regional_feature_inds{1}; vertcat( regional_feature_inds{2:4} ) }; 

% get indices for p2v values for the wfs
p2v_inds = cellfun( @(x) find(feature_labels, 'peak_to_valley', x), feature_subset_inds, 'un', 0 );
use_pca = true;

celltype_labels = eisg.celltype_class.do_shape_based_wf_class_by_subset(...
  unit_wfs, spike_labels, subset_inds, unit_wf_features, p2v_inds, n_cluster_range, use_pca);

%% Classify celltypes based on peak-to-valley across PFC and separately for BLA


%% Classify celltypes based on entire WF shape by regions


%% Classify celltypes based on entire WF across PFC and separately for BLA



%%
% 
% spike_labels = apply_cell_type_labels(...
%   spike_labels, ct_labels.label_mat, ct_labels.label_mat_cols );

function feature_string = get_feature_string_from_labels(feature_list)

feature_string = feature_list{1};
if numel( feature_list ) > 1
  for i=2:numel( feature_list )
    feature_string = [feature_string ' | ' feature_list{i}];
  end
end
feature_string = strrep(feature_string, '_', '-');

end

function pca_str = is_pca_string(use_pca)

if use_pca
  pca_str = '_pca_of_';
else
  pca_str = '';
end

end
