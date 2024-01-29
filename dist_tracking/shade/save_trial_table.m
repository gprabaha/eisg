data_dir = '/Volumes/external3/data/changlab/siqi/distance_model/Data/reformatted';

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
  
if ( 1 )
  save_p = fullfile( data_dir, 'distance_model', 'trial_table.mat' );
  shared_utils.io.require_dir( fileparts(save_p) );
  save( save_p, '-v7.3', 'trial_table' );
end
