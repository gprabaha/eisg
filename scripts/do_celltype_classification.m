%% Workspace Cleanup
clc;
clear;

 %% Loading Data
data_p = fullfile( eisg.util.project_path, 'processed_data');
fprintf('Data folder path is: %s\n', data_p);
disp('Loading data...');

% Neural data
sorted = shared_utils.io.fload( fullfile(data_p,...
  'sorted_neural_data_social_gaze.mat') );

% Behavioral data
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );

% Add task-relevant constructed ROI labels 
events = eisg.util.add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

% Celltype labels
ct_labels = shared_utils.io.fload(fullfile(data_p,...
    'celltype-labels_pfc-combined-class_p2v.mat'), 'ct_labels');
disp('Done');

%% Data Preprocessing
disp('Preprocessing data...');

% Extract 2D fcat matrix and labels from spike struct and unit WFs
[unit_spike_ts, unit_wfs, spike_labels] = eisg.util.linearize_sorted( sorted );
bfw.add_monk_labels( spike_labels ); % Add labels of participant monkeys

% Get the features of the unit's waveform (peak-te-valley, etc.)
[unit_wf_features, feature_labels] = eisg.util.get_wf_features_from_sorted( sorted );
list_of_all_features = feature_labels('wf_feature');
bfw.add_monk_labels( feature_labels );

%% Declare Parameters for the k-means Classification
disp('Declaring classification parameters...')

% Range of cluster numbers to be used for k-means classification
n_cluster_range = [2 5];

% Feature(s) to classify by
feature_list = {'peak_to_valley'};

% Do PCA on feature vectors? Function default is false
use_pca = false;

% Validity selection of units to classify
validity_filter = {'valid-unit', 'maybe-valid-unit'};

% Extract the features of interest for selected units
num_features = numel( feature_list );
features_mask = find( feature_labels, [feature_list, validity_filter] );

% Extract label indices for each region
[regional_feature_inds, regions] = findall( feature_labels, 'region', features_mask );

%% Declare Parameters for kMeans Classification
% Unit Subsets for kMeans Classification

% 1: BLA | ACC+OFC+DMPFC (used for analysis)
% 2: BLA | ACC | OFC | DMPFC
classification_subset_switch = 1;

do_classification_again = false;
do_save = false;

%% Do the Classification
switch classification_subset_switch

    case 1 % Unit Subsets for kMeans: BLA | ACC+OFC+DMPFC
        
        % Each subset specified in the cell array below is classified separately
        % This can be used to club cells of different regions together for
        % classification
        
        % Here we are classifying BLA units separately and then all PFC units
        % clubbed together separately
        subset_inds = { regional_feature_inds{1}; vertcat(regional_feature_inds{2:4}) };
        
        % Indices of the peak-to-valley values
        p2v_inds = cellfun( @(x) find(feature_labels, 'peak_to_valley', x), subset_inds, 'un', 0 );
        if do_classification_again
            disp('Starting k-means classification for subsets: BLA | ACC+OFC+DMPFC');
            [ct_labels, feature_mat]= ...
                eisg.celltype_class.do_feature_based_wf_class_for_each_subset(...
                unit_wf_features, feature_labels, subset_inds, num_features...
                , p2v_inds, n_cluster_range, use_pca);
            disp('Saving celltype classification results...');
            if do_save
                save(fullfile(data_p, 'celltype-labels_pfc-combined-class_p2v.mat'), 'ct_labels', 'feature_mat');
            end
        else
            disp('Loading saved cell-type feature matrix and labels...');
            load(fullfile(data_p, 'celltype-labels_pfc-combined-class_p2v.mat'));
        end

    case 2 % Unit Subsets for kMeans: BLA | ACC | OFC | DMPFC

        % Declare unit subsets for separate kMeans classification
        subset_inds = regional_feature_inds;

        % Indices of the peak-to-valley values
        p2v_inds = cellfun( @(x) find(feature_labels, 'peak_to_valley', x), subset_inds, 'un', 0 );
        if do_classification_again
            disp('Starting k-means classification for subsets: BLA | ACC | OFC | DMPFC');
            [ct_labels, feature_mat]= ...
                eisg.celltype_class.do_feature_based_wf_class_for_each_subset(...
                unit_wf_features, feature_labels, subset_inds, num_features...
                , p2v_inds, n_cluster_range, use_pca);
            disp('Saving celltype classification results...');
            if do_save
                save(fullfile(data_p, 'celltype-labels_region-separated-class_p2v.mat'), 'ct_labels', 'feature_mat');
            end
        else
            disp('Loading saved cell-type feature matrix and labels...');
            load(fullfile(data_p, 'celltype-labels_region-separated-class_p2v.mat'));
        end

    otherwise
        error('Unknown classification subset switch number!!');
end

%% Update Celltype Labels
disp('Updating cell-type labels...');
replace( ct_labels, 'n', 'narrow');
replace( ct_labels, 'm', 'broad');
replace( ct_labels, 'b', 'outlier');

%% Plot Celltype Waveforms for each Region

% Saving parameters
save_p = fullfile( eisg.util.project_path, 'plots/celltype_wfs', dsp3.datedir );
do_save = false;

% Add cell-type labels to a temp spike_labels variable
[uuid_I, uuid_C] = findall( ct_labels, {'uuid', 'cell-type'} ); % Gives you indices of all UUID-Celltype combinations
temp_fr_labels = spike_labels;
match_inds = bfw.find_combinations( temp_fr_labels, uuid_C(1, :) ); % Nick's bfw repo function

% Add a celltype column to temp labels
addcat( temp_fr_labels, 'cell-type' );
for i = 1:numel(match_inds)
  setcat( temp_fr_labels, 'cell-type', uuid_C{2, i}, match_inds{i} );
end

% Mask for unit validity
mask = pipe( rowmask(spike_labels) ...
  , @(m) find(spike_labels, validity_filter, m) ...
  , @(m) findnone(spike_labels, '<cell-type>', m) ...
);

% Categories within each plot
gcats = {'cell-type' };

% Categories fo different subplots
pcats = {'region'};

% If multiple figures are to be generated then find the row indices
% corresponding to each category/category combination
% We used {} since we are not separating figures by anything
fig_I = findall_or_one( spike_labels, {}, mask);
for i = 1:numel(fig_I)
  fi = fig_I{i};

  % Find the row indices of each category combination in 'pcats'
  [panel_I, panel_C] = findall( spike_labels, pcats, fi );

  % Initiate graphical objects for each unit's waveorm
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

sgtitle(['Classification by: ' get_feature_string_from_labels(feature_list)]);
if ( do_save )
  dsp3.req_savefig( gcf, save_p, prune(spike_labels(fi)), 'region', ['celltypes_classified_by' is_pca_string( use_pca ) get_feature_string_from_labels(feature_list)] );
end

%% Helper Functions

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

%% Other Classification Methods

%{
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
% do_save = true;
do_save = false;

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
sgtitle(['Classification by: ' get_feature_string_from_labels(feature_list)]);
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

ct_labels = eisg.celltype_class.do_shape_based_wf_class_by_subset(...
  unit_wfs, spike_labels, subset_inds, unit_wf_features, p2v_inds, n_cluster_range, use_pca);

%}
