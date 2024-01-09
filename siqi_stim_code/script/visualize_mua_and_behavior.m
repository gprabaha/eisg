data_root = '/Volumes/external3/data/changlab/siqi/stim';
int_p = fullfile( data_root, 'intermediates' );

spike_fs = shared_utils.io.findmat( fullfile(int_p, 'spike_ts') );
unsorted_fs = shared_utils.io.findmat( fullfile(int_p, 'unsorted_peak_voltage') );

eye_tbl = shared_utils.io.fload( fullfile(int_p, 'gaze_data_tables', 'eyes.mat') );
dot_tbl = shared_utils.io.fload( fullfile(int_p, 'gaze_data_tables', 'dot.mat') );
eye_tbl.kind = repmat( "eye", rows(eye_tbl), 1 );
dot_tbl.kind = repmat( "dot", rows(dot_tbl), 1 );
tot_tbl = [ eye_tbl; dot_tbl ];

% eg_run = "01102020_dot_4.mat";
eg_run = "01282022_position_1.mat";
sesh = string( eg_run{1}(1:8) );

raw_signal_ps = shared_utils.io.findmat( fullfile(int_p, 'pl2_mat') );
raw_signal = load( raw_signal_ps{contains(raw_signal_ps, sesh)} );

% rois
roi_ps = shared_utils.io.findmat( fullfile(int_p, 'rois') );
roi_file = shared_utils.io.fload( roi_ps{contains(roi_ps, eg_run)} );

% bounds
bounds_ps = shared_utils.io.findmat( fullfile(int_p, 'aligned_raw_samples/bounds') );
bounds = shared_utils.io.fload( bounds_ps{contains(bounds_ps, eg_run)} );

% time
time_ps = shared_utils.io.findmat( fullfile(int_p, 'aligned_raw_samples/time') );
time = shared_utils.io.fload( time_ps{contains(time_ps, eg_run)} );

% position
pos_ps = shared_utils.io.findmat( fullfile(int_p, 'aligned_raw_samples/position') );
pos = shared_utils.io.fload( pos_ps{contains(pos_ps, eg_run)} );

run_str = compose( "run_%d", str2double(eg_run{1}(max(strfind(eg_run, '_'))+1:max(strfind(eg_run, '.'))-1)) );
kind = ternary( contains(eg_run, 'dot'), "dot", "eye" );
% stim times from this session
sel = tot_tbl.kind == kind & tot_tbl.session == sesh & tot_tbl.run == run_str;

sel_tbl = tot_tbl(sel, :);

% events
evts = shared_utils.io.fload( fullfile(int_p, 'raw_events', char(eg_run)) );
vars = arrayfun( @(x) string(evts.labels(:, x)), 1:size(evts.labels, 2), 'un', 0 );
evt_labs = table( vars{:}, 'va', evts.categories );

evt_start_ts = bfw.event_column( evts, 'start_time' );
evt_stop_ts = bfw.event_column( evts, 'stop_time' );

% spikes
spikes = load( spike_fs{contains(spike_fs, sesh)} );
spike_ts = spikes.unsorted_tbl.time{1};
peak_v_file = load( unsorted_fs{contains(unsorted_fs, sesh)} );

% remove spikes within stim
stim_t0s = sel_tbl.stim_time;
ib_stim = is_ib_stim_time( spike_ts, stim_t0s, stim_t0s + 5 );
ib_thresh = peak_v_file.peak_vs < prctile( peak_v_file.peak_vs(~ib_stim), 80 );
spike_ts = spike_ts(ib_thresh);

%%

figure(1); clf; 
axs = plots.panels( [4, 1] );
hold( axs, 'on' );

% choose start / stop time, in seconds. time 0 is the start of the neural 
% recording session rather than the start of the run, so we can use
% `min(time.t)` to get the start of the run.
t0 = min( time.t ) + 60 * 1.5;  % 1.5 minutes into the run
t1 = t0 + 60;

% zoom in on specific time, if not empty
custom_xlim = []; 

% choose rois
rois = ["eyes_nf", "face"];
% rois = "eyes_nf";

%
%
% 

evt_mask = ismember(evt_labs.roi, rois) & evt_labs.looks_by == 'm1' & ...
  evt_start_ts >= t0 & evt_start_ts < t1;

y_off = 0;
evt_inds = find( evt_mask );
roi_colors = hsv( numel(rois) );
% rois
for i = 1:numel(evt_inds)
  ei = evt_inds(i);
  roi_color = roi_colors(strcmp(rois, evt_labs.roi(ei)), :);  
  ew = evt_stop_ts(ei) - evt_start_ts(ei);
  r = rectangle( 'Position', [evt_start_ts(ei), y_off, ew, y_off + 1], 'Parent', axs(1) );
  r.FaceColor = roi_color;
  r.EdgeColor = 'none';
end

% stim/sham
stim_names = ["stim", "sham"];
stim_colors = spring( numel(stim_names) );

y_off = 1;
stim_mask = stim_t0s >= t0 & stim_t0s < t1;
stim_inds = find( stim_mask );

for i = 1:numel(stim_inds)
  si = stim_inds(i);
  sw = 0.2;
  r = rectangle( 'Position', [stim_t0s(si), y_off, sw, y_off], 'Parent', axs(1) );
  r.FaceColor = stim_colors(strcmp(stim_names, sel_tbl.stim_type(si)), :);
  r.EdgeColor = 'none';
end

ylim( axs(1), [-1, y_off + 2] );
title( axs(1), 'stim + looking events' );

% signal
t0_ind = floor( t0 * raw_signal.ad.ADFreq );
t1_ind = floor( t1 * raw_signal.ad.ADFreq );
signal_t = linspace( t0, t1, t1_ind - t0_ind + 1 );
plot( axs(2), signal_t, raw_signal.ad.Values(t0_ind:t1_ind) );
title( axs(2), 'raw signal' );
ylabel( axs(2), 'voltage' );

% mua
bw_s = 0.05;
smooth_win_size_s = 300e-3;
[psth, psth_t] = bfw.event_psth( t0, {spike_ts}, {1}, {1}, 0, t1 - t0, bw_s, 'concatenate', true );

if ( 1 )
  psth = smoothdata( psth, 'gaussian', floor(smooth_win_size_s / bw_s) );
end

plot( axs(3), psth_t + t0, psth );
title( axs(3), 'mua' );
ylabel( axs(3), 'sp / s' );

% bounds
[~, bt0_ind] = min( abs(time.t - t0) );
[~, bt1_ind] = min( abs(time.t - t1) );
t_slice = time.t(bt0_ind:bt1_ind);

for i = 1:numel(rois)
  if ( 1 )
    px = pos.m1(1, bt0_ind:bt1_ind);
    py = pos.m1(2, bt0_ind:bt1_ind);
    
    curr_roi = roi_file.m1.rects(rois(i));
    expand_frac = 0.1;
    [w, h] = deal( diff(curr_roi([1, 3])), diff(curr_roi([2, 4])) );
    [x, y] = deal( mean(curr_roi([1, 3])), mean(curr_roi([2, 4])) );
    
    w = w + w * expand_frac;
    h = h + h * expand_frac;
    new_roi = [ x-w*0.5, y-h*0.5, x+w*0.5, y+h*0.5 ];
    bounds_set = shared_utils.rect.inside( new_roi, px, py );    
  else
    bounds_set = bounds.m1(rois{i});
    bounds_set = bounds_set(bt0_ind:bt1_ind);
  end
  
  [isles, durs] = shared_utils.logical.find_islands( bounds_set );
  for j = 1:numel(isles)
    x0 = t_slice(isles(j));
    w = t_slice(isles(j) + durs(j)-1) - x0;    
    hold( axs(4), 'on' );
    r = rectangle( 'Position', [x0, i, w, 1], 'Parent', axs(4) );
    r.FaceColor = roi_colors(i, :);
    r.EdgeColor = 'none';
  end
  
end
title( axs(4), 'bounds' );
ylim( axs(4), [0, numel(rois)+2] );

shared_utils.plot.match_xlims( axs );

if ( ~isempty(custom_xlim) )
  xlim( axs, custom_xlim );
end

% stim/sham text
lims = get( axs(1), 'xlim' );
min_lim = min( lims );
text( axs(1), min_lim + diff(lims) * 0.0125, 1.5, 'stim/sham' );
text( axs(1), min_lim + diff(lims) * 0.0125, 0.5, 'looking events' );

% roi text
for i = 1:numel(rois)
  lims = get( axs(4), 'xlim' );
  min_lim = min( lims );
  text( axs(4), min_lim + diff(lims) * 0.0125, i + 0.5, strrep(rois(i), '_', ' ') );  
end

%%

function ib_stim_t = is_ib_stim_time(base_ts, stim_t0s, stim_t1s)

ib_stim_t = false( numel(base_ts), 1 );
for j = 1:numel(base_ts)
  t = base_ts(j);
  if ( any(t >= stim_t0s & t <= stim_t1s) )
    ib_stim_t(j) = true;
  end
end

end