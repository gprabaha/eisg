%%

int_p = '/Volumes/external3/data/changlab/siqi/stim/intermediates';
pl2_mat_p = fullfile( int_p, 'pl2_mat' );
pl2_mats = shared_utils.io.findmat( pl2_mat_p );

% load stim/sham times
gaze_tbl = shared_utils.io.fload( fullfile(int_p, 'gaze_data_tables/eyes.mat') );

% '03172022' %  rating 3 day with 0 mua
% '07312019' %  rating 3 day with outliers
% '05272022' %  rating 3 day with good signal

pl2_p = pl2_mats{contains(pl2_mats, '02032022')};
pl2_f = shared_utils.io.filenames( pl2_p );
pl2_session = pl2_f(1:8);

match_day = string( gaze_tbl.session ) == pl2_session;
assert( sum(match_day) > 0 );

sham_ts = gaze_tbl{match_day & gaze_tbl.stim_type == 'sham', 'stim_time'};
stim_ts = gaze_tbl{match_day & gaze_tbl.stim_type == 'stim', 'stim_time'};
pl2_ad = shared_utils.io.fload( pl2_p );
filt = mua_filter( double(pl2_ad.Values), pl2_ad.ADFreq );

%%

target_intervals = [ sham_ts(:)+2, sham_ts(:)+4 ];
% target_intervals = [ stim_ts(:)+0, stim_ts(:)+4 ];
[is_mua, mua_sd] = detect_mua( filt, pl2_ad.ADFreq, target_intervals, false, 3 );

%%

meds = nan( size(target_intervals, 1), 1 );
traces = cell( size(meds) );
for j = 1:size(target_intervals, 1)
  i0 = floor( target_intervals(j, 1) * pl2_ad.ADFreq );
  i1 = floor( target_intervals(j, 2) * pl2_ad.ADFreq );
  meds(j) = median( abs(filt(i0:i1)) );
%   meds(j) = mean( abs(filt(i0:i1)) );
  traces{j, 1} = filt(i0:i1);
end
tot_trace = horzcat( traces{:} );

figure(1); clf;
subplot( 1, 2, 1 );
hist( meds, 1e3 );
title( strrep(pl2_f, '_', ' ') );
% xlim( [1e-3, 2e-3] );

subplot( 1, 2, 2 );
i0 = 0;
i0s = nan( numel(traces), 1 );
for i = 1:numel(traces)
  i0s(i) = i0;
  plot( i0:i0+numel(traces{i})-1, traces{i}, 'b' ); 
  hold on;
  i0 = i0 + numel(traces{i});
end

threshs = mua_sd * 3;
shared_utils.plot.add_horizontal_lines( gca, threshs );
shared_utils.plot.add_horizontal_lines( gca, prctile(abs(tot_trace), 99), 'r--' );
ylim( [-15e-2, 15e-2] );
% ylim( [-10, 10] );

% shared_utils.plot.add_vertical_lines( gca, i0s );