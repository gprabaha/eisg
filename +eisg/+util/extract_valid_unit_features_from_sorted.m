function valid_unit_features = extract_valid_unit_features_from_sorted(sorted_neural_data, varargin)

defaults = eisg.util.make_analysis_params_struct();
defaults.feature_list = {
    'peak_to_valley', ...
    'halfwidth', ...
    'peak_trough_ratio', ...
    'recovery_slope', ...
    'repolarization_slope'
};
defaults.use_calculated_features = true;

params = shared_utils.general.parsestruct( defaults, varargin );
feature_list = params.feature_list;
use_calculated_features = params.use_calculated_features;

valid_unit_features = struct();

feature_matrix = [];
label_matrix = [];
for file_ind = 1:size(sorted_neural_data, 2)
    n_units = sorted_neural_data(file_ind).n_units;
    for unit_ind = 1:n_units
        unit_features = [];
        if eisg.util.is_valid(sorted_neural_data, file_ind, unit_ind, params)
            for feature_ind = 1:numel(feature_list)
                if use_calculated_features
                    feature_val = sorted_neural_data(file_ind).calculated_template_features.(feature_list{feature_ind}){unit_ind};
                else
                    feature_val = sorted_neural_data(file_ind).template_features_ms4.(feature_list{feature_ind}){unit_ind};
                end
                unit_features = [unit_features feature_val];
            end
            feature_matrix = [feature_matrix; unit_features];
            uuid = sorted_neural_data(file_ind).uuid(unit_ind);
            region = sorted_neural_data(file_ind).region;
            label_matrix = [label_matrix; categorical([cellstr(region), num2str(uuid), num2str(file_ind), num2str(unit_ind)])];
        end

    end
end
label_categories = {'region', 'uuid', 'file_ind', 'unit_ind'};
valid_unit_features.feature_matrix = feature_matrix;
valid_unit_features.feature_list = feature_list;
valid_unit_features.label_matrix = label_matrix;
valid_unit_features.label_categories = label_categories;

end