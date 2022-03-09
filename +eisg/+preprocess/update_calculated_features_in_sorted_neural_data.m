function sorted_neural_data = update_calculated_features_in_sorted_neural_data(...
    sorted_neural_data, varargin)

defaults = eidg.util.make_preprocess_params_struct();
defaults.feature_list = {
    'peak_to_valley', ...
    'halfwidth', ...
    'peak_trough_ratio', ...
    'recovery_slope', ...
    'repolarization_slope'
};
defaults.n_extra_pts_per_gap_for_interpolation = 2;
defaults.include_unsure_units = true;

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
    'peak_to_valley', ...
    'halfwidth', ...
    'peak_trough_ratio', ...
    'recovery_slope', ...
    'repolarization_slope'
};

for file_ind =  1:length(sorted_neural_data)
    single_file_data = sorted_neural_data(file_ind);
    clc;
    disp('Features being calculated:');
    disp(params.feature_list);
    disp(['Starting to calculate features for file: ', single_file_data.filename]);
    disp('Feature calculation progress:');
    eidg.util.draw_progress_bar(file_ind-1, length(sorted_neural_data), params.num_ticks_in_progress_bar);
    n_units = length(single_file_data.maxchn);
    calculated_template_features_per_file = initiate_template_feature_struct(list_of_all_features, n_units);
    valid_unit_template = [];
    for unit_ind = 1:n_units
        is_unit_valid = eidg.util.is_valid(sorted_neural_data, file_ind, unit_ind, params);
        if is_unit_valid
            unit_template = single_file_data.templates(unit_ind, :);
            unit_template_norm = unit_template / max( abs(unit_template) );
            sampling_frequency = single_file_data.sort_params.sampling_frequency;
            rel_timepoints = make_timepoints_from_template(unit_template_norm, sampling_frequency);
            [interp_unit_template_norm, interp_timepts]  = interpolate_template(unit_template_norm, rel_timepoints, n_extra_pts_per_gap_for_interpolation);
            valid_unit_template(unit_ind,:) = interp_unit_template_norm;
            for feature_ind = 1:length(features_to_calculate)
                % This is where we calculate each features of a
                % particular template
                feature_name = features_to_calculate{feature_ind};
                if any( strcmp(feature_name, list_of_all_features) )
                    calculated_feature = eidg.util.calculate_single_template_feature( interp_unit_template_norm, interp_timepts, feature_name, list_of_all_features );
                    calculated_template_features_per_file.(feature_name){unit_ind} = calculated_feature;
                    
                end
            end
        end
    end
    sorted_neural_data(file_ind).calculated_template_features = calculated_template_features_per_file;
    sorted_neural_data(file_ind).normalized_templates = valid_unit_template;
    sorted_neural_data(file_ind).template_timepoints = interp_timepts;
end
clc;
disp('Features being calculated:');
disp(params.feature_list);
disp(['Starting to calculate features for file: ', single_file_data.filename]);
disp('Feature calculation progress:');
eidg.util.draw_progress_bar(file_ind, length(sorted_neural_data), params.num_ticks_in_progress_bar);

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
interp_template = interp1(timepoints, template', interp_timepts)';

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
% Scale by 1e6 and samp freq. to convert to micro s
rel_timepoints = 1e6 * rel_timepoints / sampling_frequency;

end