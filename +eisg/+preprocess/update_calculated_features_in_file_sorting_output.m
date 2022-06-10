function file_sorting_output = update_calculated_features_in_file_sorting_output(...
    file_sorting_output, varargin)

defaults = eisg.util.make_preprocess_params_struct();
defaults.feature_list = {
  'valley_ind', ...
  'peak_ind', ...
  'peak_to_valley', ...
  'halfwidth', ...
  'peak_trough_ratio', ...
  'recovery_slope', ...
  'repolarization_slope'
};
defaults.n_extra_pts_per_gap_for_interpolation = 2;
defaults.include_unsure_units = false;

params = shared_utils.general.parsestruct( defaults, varargin );

% List of template features:
% 1. peak_to_valley
% 2. halfwidth
% 3. peak_trough_ratio
% 4. recovery_slope
% 5. repolarization_slope

features_to_calculate = params.feature_list;
n_extra_pts_per_gap_for_interpolation = params.n_extra_pts_per_gap_for_interpolation;

list_of_all_features = {
  'valley_ind', ...
  'peak_ind', ...
  'peak_to_valley', ...
  'halfwidth', ...
  'peak_trough_ratio', ...
  'recovery_slope', ...
  'repolarization_slope'
};

n_units = length(file_sorting_output.maxchn);
calculated_template_features = initiate_template_feature_struct(list_of_all_features, n_units);
normalized_templates = [];
for unit_ind = 1:n_units
    unit_template = file_sorting_output.templates(unit_ind, :);
    unit_template_norm = unit_template / max( abs(unit_template) );
    sampling_frequency = file_sorting_output.sort_params.sampling_frequency;
    rel_timepoints = make_timepoints_from_template(unit_template_norm, sampling_frequency);
    [interp_unit_template_norm, interp_template_timepts]  = interpolate_template(unit_template_norm, rel_timepoints, n_extra_pts_per_gap_for_interpolation);
    normalized_templates(unit_ind,:) = interp_unit_template_norm;
    for feature_ind = 1:length(features_to_calculate)
        % This is where we calculate each features of a
        % particular template
        feature_name = features_to_calculate{feature_ind};
        if any( strcmp(feature_name, list_of_all_features) )
            calculated_feature = eisi.util.calculate_single_template_feature( interp_unit_template_norm, interp_template_timepts, feature_name, list_of_all_features );
            calculated_template_features.(feature_name){unit_ind} = calculated_feature;
        end
    end
end

file_sorting_output.calculated_template_features = calculated_template_features;
file_sorting_output.normalized_templates = normalized_templates;
file_sorting_output.template_timepoints = interp_template_timepts;

end


function calculated_template_features_per_file = initiate_template_feature_struct(list_of_all_features, n_units)

calculated_template_features_per_file = struct();
for i=1:length(list_of_all_features)
    empty_cells = cell(1, n_units);
    nan_cells = cellfun( @(x) {nan}, empty_cells );
    calculated_template_features_per_file.(list_of_all_features{i}) = nan_cells;
end

end

function [interp_template, interp_timepts]  = interpolate_template(template, timepoints, n_extra_pts_per_gap)

n_current_timepts = numel(timepoints);
n_interp_timepts = (n_current_timepts-1)*(n_extra_pts_per_gap+1)+1;
interp_timepts = linspace( timepoints(1), timepoints(end), n_interp_timepts );
interp_template = interp1(timepoints, template', interp_timepts, 'spline')';

end

function rel_timepoints = make_timepoints_from_template(unit_template, sampling_frequency)

n_points = numel( unit_template );
if rem(n_points, 2) == 1
    end_pt = (n_points - 1)/2;
    rel_timepoints = -end_pt : end_pt;
else
    end_pt = n_points / 2;
    rel_timepoints = -end_pt+1 : end_pt;
end
% Scale by samp freq. to convert to s
rel_timepoints = rel_timepoints / sampling_frequency;

end