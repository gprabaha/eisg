
%%

addpath( '/gpfs/milgram/project/chang/pg496/repositories/mpaths' );

repadd( 'chronux', true );
repadd( 'dsp3', true );
repadd( 'bfw', true );
repadd( 'eisg', true );
repadd( 'categorical', true );
repadd( 'shared_utils', true);

% if ( isempty(gcp('nocreate')) )
%   parpool( feature('numcores') );
% end

%%

data_p = '/gpfs/milgram/project/chang/pg496/';

sorted = shared_utils.io.fload( fullfile(data_p, 'sorted_neural_data_social_gaze.mat') );
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );

%%
[unit_spike_ts, unit_wfs, spike_labels] = linearize_sorted( sorted );
bfw.add_monk_labels( spike_labels );

%%
events = add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

%%

conf = bfw.set_dataroot( '/gpfs/milgram/project/chang/CHANG_LAB/naf3/Data/brains/free_viewing' );

%%
% [unit_spike_ts, unit_wfs, spike_labels] = linearize_sorted( sorted );

%%

lfp_p = bfw.gid( 'lfp', conf );

% change here to update output path
base_dst_p = '/gpfs/milgram/project/chang/pg496/social_gaze_sfcoherence_all_pairs';

lfp_files = shared_utils.io.findmat( lfp_p );
[ps, exists] = bfw.match_files( lfp_files ...
  , bfw.gid('meta', conf) ...
);


to_process = ps(all(exists, 2), :);

to_process = to_process(to_process_inds, :);

rois = { 'eyes_nf', 'face', 'right_nonsocial_object', 'right_nonsocial_object_eyes_nf_matched' };

no_nans = @(C) ~squeeze( any(any(isnan(C), 3), 2) );
not_all_nans = @(C) ~squeeze( all(all(isnan(C), 3), 2) );

%%

for i = 1:size(to_process, 1)
  fprintf( '\n %d of %d', i, size(to_process, 1) );
  
  ps = to_process(i, :);
  lfp_file = bfw.load_linked( ps{1} );
  meta_file = shared_utils.io.fload( ps{2} );
  
  %%
  
  spike_ind = find( spike_labels, meta_file.session );
  session_spike_labels = spike_labels(spike_ind);
  session_spikes = unit_spike_ts(spike_ind);
  
  lfp_regs = bfw.standardize_regions( lfp_file.key(:, 2) );
  spk_regs = bfw.standardize_regions( cellstr(session_spike_labels, 'region') );
  
  ref_ind = strcmp( lfp_regs, 'ref' );
  lfp = bfw.lfp_preprocess( lfp_file.data, 'ref_index', find(ref_ind) );
  lfp_regs = lfp_regs(~ref_ind, :);
  spks = session_spikes;
  
  lfp_labels = lfp_file.key(~ref_ind, :);
  
  %%
  
  % pairs = bfw.matching_pairs( spk_regs, lfp_regs )';
  pairs = dsp3.numel_combvec( spk_regs, lfp_regs )';
%   pairs = pairs(1:8, :);
%   pairs = at_most_n_lfp_channels( pairs, 1 );
%   pairs = pairs(1:2, :);
  
  % chan_inds = to_within_region_channel_indices( lfp_file.key );
  chan_inds = ( 1:size( lfp_file.key, 1 ) )';
  chan_inds = chan_inds(~ref_ind, :);
  lfp_labels(:, end+1) = arrayfun( @(x) sprintf('maxchn-%d', x), chan_inds, 'un', 0 );
  
  %%
  
  event_labs = events.labels';
  event_ts = bfw.event_column( events, 'start_time' );
  
  for j = 1:numel(rois)
    fprintf( '\n\t %d of %d', j, numel(rois) );
    
    event_mask = find( event_labs, [meta_file.unified_filename, rois(j)] );
    if ( isempty(event_mask) )
      continue
    end
    
    %%
    subset_event_ts = event_ts(event_mask);
    [coh, phi, f, t, info] = bfw.sfcoherence( spks, lfp, subset_event_ts, pairs ...
      , 'f_lims', [0, 85] ...
      , 'keep_if', not_all_nans ...
      , 'single_precision', true ...
    );
    
    labs = make_labels( ...
      session_spike_labels, lfp_labels, bfw.struct2fcat(meta_file), event_labs(event_mask), pairs, info.inds );
    coh = vertcat( coh{:} );    
    phi = vertcat( phi{:} );
    assert_ispair( coh, labs );
    dst_file = make_file( coh, phi, labs, f, t, info, meta_file.unified_filename );
    
    %%
    if ( 1 )
      dst_file_path = fullfile( base_dst_p, rois{j}, meta_file.unified_filename );
      shared_utils.io.require_dir( fileparts(dst_file_path) );
      shared_utils.io.psave( dst_file_path, dst_file, '-v7.3' );
    end
  end
  
end

function dst_file = make_file(coh, phi, labs, f, t, info, unified_filename)

dst_file = struct();
dst_file.unified_filename = unified_filename;
[dst_file.labels, dst_file.categories] = categorical( labs );
dst_file.coh = coh;
dst_file.phase = phi;
dst_file.f = f;
dst_file.t = t;
dst_file.info = info;

end

function l = make_labels(spk_labels, lfp_labels, meta_labels, event_labels, pairs, inds)

assert( size(pairs, 1) == numel(inds) );

l = fcat();

for i = 1:size(pairs, 1)
  spk_labs = spk_labels(pairs(i, 1));
  
  spk_chan = cellstr( spk_labs, 'maxchn' );
  lfp_chan = lfp_labels(pairs(i, 2), 1);
  spk_lfp_chan = strjoin( [spk_chan, lfp_chan], '_' );
  addsetcat( spk_labs, 'channel', spk_lfp_chan );
  
  spk_reg = cellstr( spk_labs, 'region' );
  lfp_reg = lfp_labels(pairs(i, 2), 2);
  lfp_max_chn = lfp_labels(pairs(i, 2), 3);
  spk_lfp_reg = strjoin( [spk_reg, lfp_reg], '_' );
  setcat( spk_labs, 'region', spk_lfp_reg );
  addsetcat( spk_labs, 'lfp-maxchn', sprintf('lfp-%s', char(lfp_max_chn)) );
  
  join( spk_labs, meta_labels );
  li = join( event_labels', spk_labs );  
  li = li(inds{i});  
  
  if ( i == 1 )
    l = li;
  else
    append( l, li );
  end
end

end

function pairs = at_most_n_lfp_channels(pairs, n)

[~, ~, ic] = unique( pairs(:, 1) );
ic = groupi( ic );
ic = cellfun( @(x) x(1:min(n, numel(x))), ic, 'un', 0 );
pairs = pairs(vertcat(ic{:}), :);

end

function inds = to_within_region_channel_indices(chan_matrix)

[regs, ~, ic] = unique( chan_matrix(:, 2) );
ic = groupi( ic );
inds = zeros( size(chan_matrix, 1), 1 );
for i = 1:numel(ic)
  inds(ic{i}) = 1:numel(ic{i});
end

end