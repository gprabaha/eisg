function [mdls, cell_labels, fit_func_name] = run_distance_model(spike_file, trial_table, all_spike_counts_in_fixation_interval, use_discretized)

num_cells = numel( all_spike_counts_in_fixation_interval );
cell_labels = nan( num_cells, 1 );
mdls = cell( num_cells, 1 );

% fit_fn = @fit_step_wise_distance_model;
fit_fn = @fit_non_step_wise_distance_model;

fit_func_name = func2str( fit_fn );

distances_to_m2s_eyes = trial_table.distances_to_m2s_eyes;
m2_dist_to_m1s_eyes = trial_table.m2_dist_to_m1s_eyes;
m2_dist_is_valid = trial_table.m2_dist_is_valid;
m1_fix_durations = trial_table.m1_fix_durations;
sessions = trial_table.sessions;
contra_ipsi_mask = trial_table.contra_ipsi_mask;

parfor cell_index = 1:num_cells

fprintf( '\n %d of %d', cell_index, num_cells );

spike_session = char( spike_file.labels(cell_index, 'session') );
index_of_fixations_this_session = sessions == spike_session;
index_of_fixations_this_session = index_of_fixations_this_session & contra_ipsi_mask;

self_distance = distances_to_m2s_eyes(index_of_fixations_this_session);
self_distance(self_distance > 20) = nan;

other_distance = m2_dist_to_m1s_eyes(index_of_fixations_this_session);
other_is_valid = m2_dist_is_valid(index_of_fixations_this_session);
other_distance(~other_is_valid) = nan;

if ( use_discretized )
    bin_edges = [0, 5, 10, 15, 20.0001];
    self_distance = discretize( self_distance, bin_edges );
    other_distance = discretize( other_distance, bin_edges );
end

y = all_spike_counts_in_fixation_interval{cell_index};
%   double check whether units of 'offset' should be ms or s, they are
%   currently s
offset = m1_fix_durations(index_of_fixations_this_session) / 1e3;
X = [self_distance, other_distance];

only_fit = ~any( isnan(X), 2 ) & ~isnan( y ) & ~isnan( offset );
mdl = fit_fn( X(only_fit, :), y(only_fit), offset(only_fit) );

cell_labels(cell_index) = determine_cell_model_label( mdl );
mdls{cell_index} = mdl;

%cell_labels(1:10)

end

end

function mdl = fit_non_step_wise_distance_model(X, y, offset)

mdl = fitglm(...
    X, y, 'linear', 'Link','log','Distribution' ,'poisson','Offset',log(offset) );

end

function mdl = fit_step_wise_distance_model(X, y, offset)

mdl = stepwiseglm(...
    X, y, 'constant','Lower','linear','Upper','linear' ...
    , 'Link','log','Distribution' ,'poisson','Offset',log(offset), 'Verbose', 0);

end

%   from the model output, label the cell with one of 4 possible labels.
%   label 1 means only the x1 term was significant in the model
%       x1 is the term corresponding to self distance
%   label 2 means only the x2 term was significant in the model
%       x2 is the term corresponding to other distance
%   label 3 means both terms were significant
%   label 4 means neither term was significant
function label = determine_cell_model_label(mdl)



x1_sig = is_term_significant( mdl, 'x1' );
x2_sig = is_term_significant( mdl, 'x2' );

if (x1_sig > 0 & x2_sig == 0);
    label = 1; % label should be 1, 2, 3 or 4
elseif x1_sig == 0 & x2_sig > 0;
    label = 2;
elseif x1_sig > 0 & x2_sig > 0;
    label = 3; 
elseif x1_sig == 0 & x2_sig == 0;
    label = 4;

d = 10;

end

d= 10;
end

%   to determine whether a term is significant:
%   1. first check to see whether the term is in the coefficients table
%       in the model. If it isn't, then the term is not significant.
%   2. if it is in the coefficients table, then check to see if the
%       associated pValue is less than 0.05 . If not, then the term is not
%       significant.
%   3. otherwise, if the term is in the table and has a pValue < 0.05, the
%       term is significant
function tf = is_term_significant(mdl, term)

tf = false;

coefficients = mdl.Coefficients;
coefficient_names = coefficients.Properties.RowNames;
coefficient_pvalues = mdl.Coefficients.pValue;

[is_in_set, loc] = ismember(term, coefficients.Properties.RowNames);
if ~is_in_set 
    tf = false;
    return
end

p_value_for_term = mdl.Coefficients.pValue(loc);
tf = p_value_for_term < 0.05;




x = 10;

%   hint: use ismember to check whether term is in set of coefficient names
%   is_in_set = ...
%   if ~is_in_set % term is not in model, so can't be significant
%       tf = false;
%       return
%   end
%
%   p_value_for_term = ...
%   tf = p_value_for_term < 0.05

end