clc;
clear;

data_p = fullfile( eisg.util.project_path, 'processed_data');

sorted = shared_utils.io.fload( fullfile(data_p, 'sorted_neural_data_social_gaze.mat') );
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );
ct_labels = shared_utils.io.fload(fullfile(data_p, 'celltype-labels_pfc-combined-class_p2v.mat'), 'ct_labels'); % Celltype labels ct_labels

%%

[unit_spike_ts, unit_wfs, spike_labels] = eisg.util.linearize_sorted( sorted );
bfw.add_monk_labels( spike_labels );

[uuid_I, uuids] = findall( spike_labels, 'uuid', find(spike_labels, {'valid-unit', 'maybe-valid-unit'}) );
match_I = bfw.find_combinations( ct_labels, uuids );

for i = 1:numel(uuid_I)
  if ( ~isempty(match_I{i}) )
    ct_label = cellstr( ct_labels, 'cell-type', match_I{i} );
    addsetcat( spike_labels, 'cell-type', ct_label, uuid_I{i} );
  end
end

events = eisg.util.add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

%%
rois = { 'whole_face', 'right_nonsocial_object_whole_face_matched', 'eyes_nf', 'face' };

evt_mask = find( events.labels, [{'m1'}, rois] );
spk_mask = find( spike_labels, {'valid-unit', 'maybe-valid-unit'} );

min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;

[psth_matrix, psth_labels, t] = compute_psth(...
    unit_spike_ts, spike_labels, spk_mask ...
  , evts, events.labels, evt_mask ...
  , min_t, max_t, bin_width );

psth_matrix = smoothdata( psth_matrix, 2, 'movmean', 10 );

%% Ranksum comparison for all bins of all unit

num_bins = size( psth_matrix, 2 );
% [~,uuids] = findall(psth_labels, 'uuid');
% ps = nan( numel( uuids ), 1 );
% ps_e_nef = nan( numel( uuids ), 1 );
ps_social_nonsocial = nan( numel( uuids ), num_bins );
ps_features = nan( numel( uuids ), num_bins );
parfor i=1:num_bins
%   rs_outs{i} = dsp3.ranksum( psth_matrix(:, i), psth_labels, {'uuid'}, 'whole_face', 'right_nonsocial_object_whole_face_matched' );
%   rs_outs_e_nef{i} = dsp3.ranksum( psth_matrix(:, i), psth_labels, {'uuid'}, 'face', 'eyes_nf' );
  
  rs_outs_social_nonsocial{i} = dsp3.ranksum( psth_matrix(:, i), psth_labels, {'uuid'}, 'whole_face', 'right_nonsocial_object_whole_face_matched' );
  rs_outs_feature{i} = dsp3.ranksum( psth_matrix(:, i), psth_labels, {'uuid'}, 'eyes_nf', 'face' );
  
%   ps(:,i) = cellfun( @(x) x.p, rs_outs{i}.rs_tables );
%   ps_e_nef(:,i) = cellfun( @(x) x.p, rs_outs_e_nef{i}.rs_tables );
  ps_social_nonsocial(:,i) = cellfun( @(x) x.p, rs_outs_social_nonsocial{i}.rs_tables );
  ps_features(:,i) = cellfun( @(x) x.p, rs_outs_feature{i}.rs_tables );
  
  fprintf( 'compared bin %d of %d\n', i, num_bins );
end

%%

p_thresh = 0.05;
seq_dur_thresh = 5;

sig_social_color = [1, 0, 0]; % red
sig_eye_color = [0, 1, 0]; % green
both_color = [0, 0, 1]; % blue

image_array = ones( size(ps_social_nonsocial, 1), size(ps_social_nonsocial, 2), 3 );

sig_social = ps_social_nonsocial < p_thresh;
sig_eye = ps_features < p_thresh;
sig_both =  sig_eye & sig_social;


%image_array = assign_at( image_array, sig_social & ~sig_eye, sig_social_color );
image_array = assign_at( image_array, sig_social , sig_social_color );
%image_array = assign_at( image_array, sig_eye & ~sig_social, sig_eye_color );
%image_array = assign_at( image_array, sig_eye & sig_social, both_color );

% tot_sig_social_bins = sum(sig_social,2);
% tot_sig_eye_bins = sum(sig_eye,2);
% tot_bins = tot_sig_social_bins + tot_sig_eye_bins;
% [~,sorted_inds] = sort(tot_bins);
% image_array_sorted = image_array(sorted_inds, :, :);

%unit_mask = findnone( rs_outs_social_nonsocial{1}.rs_labels, {'b'} );
unit_mask = findnone( rs_outs_social_nonsocial{1}.rs_labels, {'b', 'ofc', 'dmpfc'} );

[I, C] = findall( rs_outs_social_nonsocial{1}.rs_labels, {'cell-type', 'region'}, unit_mask );
axs = plots.panels( numel(I) );
for i = 1:numel(axs)
  %ind = arrayfun( @(x)( find(sorted_inds==x) ), I{i} );
  ind = I{i};
  % Count of ind is the number of cells of that type in that region
  ind2 = [];
  for j = ind'
    [~, seq_durs] = shared_utils.logical.find_islands( sig_social(j,:) );
    if max(seq_durs)>=seq_dur_thresh
      ind2 = [ind2; j];
    end
  end
  ind = ind2;
  tot_sig_social_bins = sum(sig_social(ind,:),2);
  tot_sig_eye_bins = sum(sig_eye(ind,:),2);
  tot_sig_both = sum(sig_both(ind,:),2);
  tot_bins = tot_sig_social_bins + tot_sig_eye_bins;
  %[~,sorted_ind] = sort(tot_bins);
  [~,sorted_ind] = sort(tot_sig_social_bins);
  %[~,sorted_ind] = sort(tot_sig_eye_bins);
  %[~,sorted_ind] = sort(tot_sig_both);
  image_array_sorted = image_array(ind(sorted_ind), :, :);
  plot( axs(i), t, sum(sig_social(ind,:)) );
  t_pre = t>-0.5 & t<=0;
  t_post = t>0 & t<=0.5;
  tot_cells_per_bin_pre = sum(sig_social(ind,t_pre));
  tot_cells_per_bin_post = sum(sig_social(ind,t_post));
  median_pre_cells = median(tot_cells_per_bin_pre);
  median_post_cells = median(tot_cells_per_bin_post);
  mean_pre_cells = mean(tot_cells_per_bin_pre);
  mean_post_cells = mean(tot_cells_per_bin_post);
  ranksum_p = ranksum(tot_cells_per_bin_pre, tot_cells_per_bin_post);
  [~, ttest_p] = ttest2(tot_cells_per_bin_pre, tot_cells_per_bin_post);
  celltype_and_region = char(fcat.strjoin(C(:, i), ' | '));
  fprintf('%s | median #cells/bin pre: %0.1f; post: %0.1f; ranksum p: %0.3f;\n', celltype_and_region, median_pre_cells, median_post_cells, ranksum_p);
  fprintf('%s | mean #cells/bin pre: %0.2f; post: %0.2f; 2 sample t-test p: %0.3f;\n', celltype_and_region, mean_pre_cells, mean_post_cells, ttest_p);
  %histogram( axs(i), sum(sig_social(ind,:), 2) );
  %plot( axs(i), t, sum(sig_eye(ind,:)) );
  %imagesc( axs(i), t, 1:numel(ind), image_array_sorted );
  title( axs(i), strrep(fcat.strjoin(C(:, i), ' | '), '_', ' ') );
end

% Try isolating islands of at least 5 or so bins, and then sorting by the location of the
% first island

% imagesc( t, 1:size(image_array, 2), image_array );

%% AUC analysis


%[auc_labels, unit_I] = keepeach( psth_labels', 'uuid', find(psth_labels, ref(combs(psth_labels, 'uuid'), '()', 1)) );
[auc_labels, unit_I] = keepeach( psth_labels', 'uuid' );
roi_a = 'whole_face';
roi_b = 'right_nonsocial_object_whole_face_matched';

aucs = nan( numel(unit_I), size(psth_matrix, 2) );
null_aucs = cell( numel(unit_I), 1 );
perm_p = nan( size(aucs) );
z_scored_aucs = nan( size(aucs) );

for i = 1:numel(unit_I)
  fprintf( '%d of %d\n', i, numel(unit_I) );
  ind_a = find( psth_labels, roi_a, unit_I{i} );
  ind_b = find( psth_labels, roi_b, unit_I{i} );
  aucs(i, :) = auc_over_time( psth_matrix, ind_a, ind_b ); 
  null_aucs{i} = auc_perm_test( psth_matrix, ind_a, ind_b, 100 );
  
  z_scored_aucs(i, :) = (aucs(i, :) - (mean(null_aucs{i}, 1))) ./ std( null_aucs{i}, [], 1 );
  
  %perm_p(i, :) = sum( aucs(i, :) < null_aucs{i}, 1 ) ./ 100;
end

%%
unit_mask = findnone( auc_labels, {'b', 'ofc', 'dmpfc'} );

[I, C] = findall( auc_labels, {'cell-type', 'region'}, unit_mask );
axs = plots.panels( numel(I) );

% Sort the array by the location of max deviation from 0.5

for i = 1:numel(axs)
  %ind = arrayfun( @(x)( find(sorted_inds==x) ), I{i} );
  ind = I{i};
  tot_sig_social_bins = sum(sig_social(ind,:),2);
  [~,sorted_ind] = sort(tot_sig_social_bins);
  %imagesc( axs(i), t, 1:numel(ind), aucs(ind(sorted_ind), :) );
  imagesc( axs(i), t, 1:numel(ind), z_scored_aucs(ind(sorted_ind), :) );
  colorbar( axs(i) );
  set( axs(i), 'clim', [-8, 8] );
  %set( axs(i), 'clim', [0.3, 0.65] );
  title( axs(i), strrep(fcat.strjoin(C(:, i), ' | '), '_', ' ') );
end
%

%%
sig_cell_bins_e = zeros( size( ps_social_nonsocial ) );
sig_cell_bins_e(ps_social_nonsocial < 0.05) = 1;

sig_cell_bins_f = zeros( size( ps_features ) );
sig_cell_bins_f(ps_features < 0.05) = 1;

%%

islands = {};
for i = 1:size(sig_social, 1)
  [seq_starts, seq_durs] = shared_utils.logical.find_islands( sig_social(i,:) );
  islands{i} = [seq_starts', seq_durs];
%   if ( ~isempty(seq_durs) )
%     max_consec_bins(i) = max( seq_durs );
%   end
%   sig_bins = find( sig_cell_bins(i,:) );
%   if ~isempty( sig_bins )
%     consec_bins = 1;
%     max_consec_bins(i) = 1;
%     sig_bin_diff = diff( sig_bins );
%     for j = 1:numel(sig_bin_diff)
%       if sig_bin_diff(j) == 1
%         consec_bins = consec_bins + 1;
%         if consec_bins > max_consec_bins(i)
%           max_consec_bins(i) = consec_bins;
%         end
%       else
%         consec_bins = 1;
%       end
%     end
%   end
end


%%

sig_cells = max_consec_bins >= 5;
use_labels = rs_outs{1}.rs_labels';
[cell_type_I, cell_types] = findall( use_labels, {'cell-type', 'region'}, findnone(use_labels, 'b') );
%cell_type_props = cellfun( @(x) sum(sig_cells(x)) / numel(sig_cells(x)), cell_type_I );
cell_type_props = cellfun( @(x) pnz(sig_cells(x)), cell_type_I );

prop_table = table( cell_type_props, 'RowNames', fcat.strjoin(cell_types, '_') );
prop_table

[tbl, reg_labels] = do_prop_tests( sig_cells, use_labels, findnone(use_labels, 'b') );
tbl


%%  plot single unit psth, overlaying bins with significant discrimination

is_sig_e = sig_cell_bins_e;
is_sig_f = sig_cell_bins_f;
assert_ispair( ps, use_labels );

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
      sig_sub = find( is_sig_e(match_sig, :) );
      plot( gca, t(sig_sub), repmat(max(get(gca, 'ylim')), numel(sig_sub), 1), 'k*' );
      sig_sub_e_nef = find( is_sig_f(match_sig, :) );
      plot( gca, t(sig_sub_e_nef), repmat(0.97*max(get(gca, 'ylim')), numel(sig_sub_e_nef), 1), 'r*' );
    end
    
    rep_join = @(s) strrep( fcat.strjoin(s, ' | '), '_', ' ' );
    
    legend( hs, rep_join(g_C), 'Location', 'Best' );
    title( rep_join(p_C(:, j)) );
  end
  
  if ( 1 )
    reg_subdir = fig_C{2, i};
    save_labs = prune( psth_labels(fi) );
    save_p = fullfile( eisg.util.project_path(), 'data/plots/psth_ranksum_diff', dsp3.datedir );
    save_p = fullfile( save_p, reg_subdir );
    dsp3.req_savefig( gcf, save_p, save_labs, 'uuid' );
  end
end


%%

function [tbl, reg_labels, tbls] = do_prop_tests(meets_criterion, labels, mask)

assert_ispair( meets_criterion, labels );

[reg_labels, reg_I] = keepeach( labels', 'region', mask );
n_I = cellfun( @(x) find(labels, 'n', x), reg_I, 'un', 0 );
m_I = cellfun( @(x) find(labels, 'm', x), reg_I, 'un', 0 );

tbls = do_prop_test( meets_criterion, n_I, m_I );
chi2ps = [ tbls.p ]';
tbl = array2table( chi2ps, 'RowNames', fcat.strjoin(reg_labels(:, 'region')') );

end

function tbls = do_prop_test(meets_criterion, inds_a, inds_b)

xa = cellfun( @(x) sum(meets_criterion(x)), inds_a );
xb = cellfun( @(x) sum(meets_criterion(x)), inds_b );
na = cellfun( @numel, inds_a );
nb = cellfun( @numel, inds_b );

tbls = struct( 'h', {}, 'p', {}, 'chi2stat', {} );
for i = 1:numel(xa)
  [h, p, chi2stat] = prop_test( [xa(i), xb(i)], [na(i), nb(i)], false );
  tbls(i) = struct( 'h', h, 'p', p, 'chi2stat', chi2stat );
end

end

function im = assign_at(im, inds, color)

[row, col] = find( inds );

for i = 1:numel(row)
  im(row(i), col(i), :) = color;
end

end


function aucs = auc_over_time(spikes, ind_a, ind_b)

aucs = nan( 1, size(spikes, 2) );

for i = 1:size(spikes, 2)
  spks_a = spikes(ind_a, i);
  spks_b = spikes(ind_b, i);
  auc = score_auc( spks_a, spks_b );
  aucs(i) = auc;
end

end

function null_aucs = auc_perm_test(spikes, ind_a, ind_b, perm_iters)

null_aucs = cell( perm_iters, 1 );
for i = 1:perm_iters
  [ia, ib] = shuffle2( ind_a, ind_b );
  null_aucs{i} = auc_over_time( spikes, ia, ib );
end
null_aucs = vertcat( null_aucs{:} );
null_aucs = sort( null_aucs );

end

function [ic, id] = shuffle2(ia, ib)

i = [ia; ib];
i = i(randperm(numel(i)));
ic = i(1:numel(ia));
id = i(numel(ia)+1:end);
assert( numel(ic) == numel(ia) && numel(id) == numel(ib) );

end

function auc = score_auc(a, b)

t = false( numel(a) + numel(b), 1 );
t(1:numel(a)) = true;
y = [ a; b ];
auc = scoreAUC( t, y );

end


%% This is the fast AUC function copied over
% https://www.mathworks.com/matlabcentral/fileexchange/50962-fast-auc

function auc = scoreAUC(labels,scores)
% Calculates the AUC - area under the curve.
%
% Besides being the area under the ROC curve, AUC is has a slightly 
% less known interpretation:
% If you choose a random pair of samples which is one positive and one
% negative - AUC is the probabilty that the positive-sample score is above
% the negative-sample score.
% Here we compute the AUC by counting these pairs.
% 
% auc = scoreAUC(labels,scores)
% N x 1 boolean labels
% N x 1 scores
%
% http://www.springerlink.com/content/nn141j42838n7u21/fulltext.pdf
%
% ==== Lior Kirsch 2014 ====

assert( islogical(labels) ,'labels input should be logical');
assert( isequal(size(labels),size(scores)) ,'labels and scores should have the same size');

[n,m] = size(labels);
assert( m==1, 'should have sizse of (n,1)');

num_pos = sum(labels);
num_neg = n - num_pos;
assert( ~(num_pos==0), 'no positive labels entered');
assert( ~(num_pos==n), 'no negative labels entered');


ranks = tiedrank(scores);
auc = ( sum( ranks(labels) ) - num_pos * (num_pos+1)/2) / ( num_pos * num_neg);

end