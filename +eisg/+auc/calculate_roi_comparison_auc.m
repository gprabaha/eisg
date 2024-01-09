function [aucs, z_scored_aucs, auc_labels] = calculate_roi_comparison_auc(...
    psth_matrix, psth_labels, roi_a, roi_b...
    )

if nargin < 4
    roi_a = 'whole_face';
    roi_b = 'right_nonsocial_object_whole_face_matched';
end
disp('Starting AUC computations...');
[auc_labels, unit_I] = keepeach( psth_labels', 'uuid' );
aucs = nan( numel(unit_I), size(psth_matrix, 2) );
z_scored_aucs = nan( size(aucs) );
fprintf('Calculating for unit: ');
parfor i = 1:numel(unit_I)
    if mod(i, 20) ~= 0
        fprintf('%d ', i);
    else
        fprintf('%d \n', i);
    end
    drawnow;
    ind_a = find( psth_labels, roi_a, unit_I{i} );
    ind_b = find( psth_labels, roi_b, unit_I{i} );
    aucs(i, :) = auc_over_time( psth_matrix, ind_a, ind_b ); 
    null_aucs = auc_perm_test( psth_matrix, ind_a, ind_b, 100 );
    z_scored_aucs(i, :) = (aucs(i, :) - mean( null_aucs, 1 )) ./ std( null_aucs, [], 1 );
end
disp('Done');

end

%% Helper Functions

function aucs = auc_over_time(psth_matrix, ind_a, ind_b)
aucs = nan( 1, size(psth_matrix, 2) );
for i = 1:size(psth_matrix, 2)
    spks_a = psth_matrix(ind_a, i);
    spks_b = psth_matrix(ind_b, i);
    auc = score_auc( spks_a, spks_b );
    aucs(i) = auc;
end
end

function auc = score_auc(a, b)
t = false( numel(a) + numel(b), 1 );
t(1:numel(a)) = true;
y = [ a; b ];
auc = scoreAUC( t, y );
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

% 2D Shuffle
function [ic, id] = shuffle2(ia, ib)
i = [ia; ib];
i = i(randperm(numel(i)));
ic = i(1:numel(ia));
id = i(numel(ia)+1:end);
assert( numel(ic) == numel(ia) && numel(id) == numel(ib) );
end


%% Downloaded Functions
% Fast AUC
% Copied over from: https://www.mathworks.com/matlabcentral/fileexchange/50962-fast-auc

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