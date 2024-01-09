
%%
rois = { 'whole_face', 'right_nonsocial_object_whole_face_matched',...
  'eyes_nf', 'face' };
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


%%

[auc_labels, unit_I] = keepeach( psth_labels', 'uuid' );
roi_a = 'whole_face';
roi_b = 'right_nonsocial_object_whole_face_matched';
aucs = nan( numel(unit_I), size(psth_matrix, 2) );
null_aucs = cell( numel(unit_I), 1 );
perm_p = nan( size(aucs) );
z_scored_aucs = nan( size(aucs) );
parfor i = 1:numel(unit_I)
  % fprintf( '%d of %d\n', i, numel(unit_I) );
  ind_a = find( psth_labels, roi_a, unit_I{i} );
  ind_b = find( psth_labels, roi_b, unit_I{i} );
  aucs(i, :) = auc_over_time( psth_matrix, ind_a, ind_b ); 
  null_aucs{i} = auc_perm_test( psth_matrix, ind_a, ind_b, 100 );
  z_scored_aucs(i, :) = (aucs(i, :) - mean( null_aucs{i}, 1 )) ./ std( null_aucs{i}, [], 1 );  
end

%%

t_snip_fix = t>-0.25 & t<=0.25;
t_snip_flank = (t>-0.5 & t<=-0.25) | (t>0.25 & t<=0.5);

t_snip_pre = t>-0.5 & t<=0;
t_snip_post = t>0 & t<=0.5;

i=1;

abs_zscore_auc = [abs_zscore_auc_broad abs_zscore_auc_narrow];

n = numel(abs_zscore_auc{i}(:,1));
for j = 1:n
  if ~isempty( rmmissing( mean( abs_zscore_auc{i}(j, t_snip_fix) ) ) ) && ~isempty( rmmissing( mean( abs_zscore_auc{i}(j, t_snip_flank) ) ) )
    p_fix_flank(j) = ranksum( abs_zscore_auc{i}(j, t_snip_fix), abs_zscore_auc{i}(j, t_snip_flank) );
  else
    p_fix_flank(j) = nan;
  end
  if ~isempty( rmmissing( abs_zscore_auc{i}(j, t_snip_pre) ) )  && ~isempty( rmmissing( mean( abs_zscore_auc{i}(j, t_snip_post) ) ) )
    p_pre_post(j) = ranksum( abs_zscore_auc{i}(j, t_snip_pre), abs_zscore_auc{i}(j, t_snip_post) );
  else
    p_pre_post(j) = nan;
  end
end

X1 = nansum( p_fix_flank<0.05 );
Y1 = nansum( p_pre_post<0.05 );
N1 = n;

clear p_fix_flank p_pre_post;

i=i+1;
n = numel(abs_zscore_auc{i}(:,1));
for j = 1:n
  if ~isempty( rmmissing( mean( abs_zscore_auc{i}(j, t_snip_fix) ) ) ) && ~isempty( rmmissing( mean( abs_zscore_auc{i}(j, t_snip_flank) ) ) )
    p_fix_flank(j) = ranksum( abs_zscore_auc{i}(j, t_snip_fix), abs_zscore_auc{i}(j, t_snip_flank) );
  else
    p_fix_flank(j) = nan;
  end
  if ~isempty( rmmissing( mean( abs_zscore_auc{i}(j, t_snip_pre) ) ) ) && ~isempty( rmmissing( mean( abs_zscore_auc{i}(j, t_snip_post) ) ) )
    p_pre_post(j) = ranksum( abs_zscore_auc{i}(j, t_snip_pre), abs_zscore_auc{i}(j, t_snip_post) );
  else
    p_pre_post(j) = nan;
  end
end

X2 = nansum( p_fix_flank<0.05 );
Y2 = nansum( p_pre_post<0.05 );
N2 = n;

disp([Y1 Y2 N1 N2]);
[~, p_flank] = prop_test([X1 X2], [N1 N2], false)
[~, p_pre_post] = prop_test([Y1 Y2], [N1 N2], false)


%%

function im_out = assign_at(im, inds, color)

im_out = im;
[row, col] = find( inds );
for i = 1:numel(row)
  im_out(row(i), col(i), :) = color;
end

end

% AUC Over Time
function aucs = auc_over_time(spikes, ind_a, ind_b)

aucs = nan( 1, size(spikes, 2) );
for i = 1:size(spikes, 2)
  spks_a = spikes(ind_a, i);
  spks_b = spikes(ind_b, i);
  auc = score_auc( spks_a, spks_b );
  aucs(i) = auc;
end

end

% AUC Perm Test
function null_aucs = auc_perm_test(spikes, ind_a, ind_b, perm_iters)

null_aucs = cell( perm_iters, 1 );
for i = 1:perm_iters
  [ia, ib] = shuffle2( ind_a, ind_b );
  null_aucs{i} = auc_over_time( spikes, ia, ib );
end
null_aucs = vertcat( null_aucs{:} );
null_aucs = sort( null_aucs );

end

% 2D Shuffle
function [ic, id] = shuffle2(ia, ib)

i = [ia; ib];
i = i(randperm(numel(i)));
ic = i(1:numel(ia));
id = i(numel(ia)+1:end);
assert( numel(ic) == numel(ia) && numel(id) == numel(ib) );

end

% AUC Score
function auc = score_auc(a, b)

t = false( numel(a) + numel(b), 1 );
t(1:numel(a)) = true;
y = [ a; b ];
auc = scoreAUC( t, y );

end

%%Fast AUC
%%Copied over from: https://www.mathworks.com/matlabcentral/fileexchange/50962-fast-auc

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
