%%  load in data

data_dir = '/Volumes/external3/data/changlab/siqi/distance_model/Data/reformatted';

trial_table = shared_utils.io.fload( fullfile(data_dir, 'distance_model/trial_table.mat') );
spike_file = shared_utils.io.fload( fullfile(data_dir, 'relabeled_cells.mat') );
bfw.add_monk_labels( spike_file.labels );
spike_labels = fcat.totable( spike_file.labels );

%%  compute psth for all cells

sessions = trial_table.sessions;
mdl_prefix = 'subsampled_ipsi';
mdl_prefix = 'contra';
mdl_prefix = 'ipsi';
% mdl_prefix = 'all_contra_ipsi'; % everything

shuffle_prefix = '';
contra_ipsi_mask = make_real_contra_ipsi_mask( trial_table.m1_hemifield_origin_delta_pos, sessions, mdl_prefix );

%%  compute psth using some subset of trials given by contra_ipsi_mask

include_psth_timecourse = false;

[allpsth, psth_ts, all_spike_counts_in_fixation_interval, index_of_fixations] = ...
    compute_contra_ipsi_psth( spike_file, sessions ...
    , contra_ipsi_mask, trial_table.fixation_start_ts ...
    , trial_table.fixation_stop_ts, include_psth_timecourse );

assert( numel(allpsth) == numel(spike_file.spike_times));

%%  run distance model for all cells

model_table = trial_table;
model_table.contra_ipsi_mask = contra_ipsi_mask;

use_discretized = false;

discretized_prefix = '';
if ( use_discretized )
  discretized_prefix = 'disc-';
end

[mdls, cell_labels, fit_function_name] = run_distance_model( ...
  spike_file, model_table, all_spike_counts_in_fixation_interval, use_discretized );

if ( 1 )
  file_name = sprintf( '%s%s%s-model.mat', discretized_prefix, mdl_prefix, shuffle_prefix );
  save_p = fullfile( data_dir, 'distance_model', dsp3.datedir, file_name );
  shared_utils.io.require_dir( fileparts(save_p) );
  save( save_p, "mdls", "cell_labels", "spike_labels", "fit_function_name" );
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
%         contra_ipsi_mask = m1_hemifield_origin_delta_pos(:, 1);
        contra_ipsi_mask(:) = true;
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