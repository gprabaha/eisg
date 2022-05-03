function [model_p, real_p, null_ps] = binary_fitcdiscr_perm_test(data, ind_a, ind_b, varargin)

defaults = struct();
defaults.hold_out = 0.1;
defaults.iters = 1e3;
params = shared_utils.general.parsestruct( defaults, varargin );

hold_out = params.hold_out;
iters = params.iters;

real_p = one_iter( data, ind_a, ind_b, hold_out );
null_ps = zeros( iters, 1 );

for i = 1:iters      
  inds = [ind_a; ind_b];
  inds = inds(randperm(numel(inds)));
  null_ind_a = inds(1:numel(ind_a));
  null_ind_b = inds(numel(ind_a)+1:end);
  null_ps(i) = one_iter( data, null_ind_a, null_ind_b, hold_out );
end

model_p = 1 - pnz( real_p > null_ps );

end

function p_match = one_iter(data, ind_a, ind_b, hold_out)

keep_in = 1 - hold_out;
num_a = max( 1, floor(numel(ind_a) * keep_in) );
num_b = max( 1, floor(numel(ind_b) * keep_in) );

ind_aa = ind_a(randperm(numel(ind_a), num_a));
ind_bb = ind_b(randperm(numel(ind_b), num_b));

[train_Y, train_I] = combine_indices( ind_aa, ind_bb );
train_X = data(train_I);

rest_a = setdiff( ind_a, ind_aa );
rest_b = setdiff( ind_b, ind_bb );
[test_Y, test_I] = combine_indices( rest_a, rest_b );
test_X = data(test_I);

res = fitcdiscr( train_X, train_Y );
pred = predict( res, test_X );
p_match = sum( pred == test_Y ) / numel( test_Y );

end

function [y, i] = combine_indices(a, b)
i = [ a; b ];
y = [ ones(numel(a), 1); zeros(numel(b), 1) ];
end