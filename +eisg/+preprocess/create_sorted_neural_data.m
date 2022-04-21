function sorted_neural_data = create_sorted_neural_data(varargin)

defaults = eisg.util.make_preprocess_params_struct();
defaults.path_to_sorter_output = '/Volumes/ExtSSD/sorted_neural_data/social_gaze/';
defaults.integrate_validity = true;
defaults.path_to_validity = '/Users/prabaha/repositories/eidg/processed_data/unit_validity_social_gaze.mat';
defaults.calculate_features = true;
defaults.include_unsure_units = true;

params = shared_utils.general.parsestruct( defaults, varargin );

validity = load( params.path_to_validity );
validity = validity.validity;

main_dir = params.path_to_sorter_output;
all_folders = dir(main_dir);
all_folders = all_folders([all_folders.isdir]);
all_folders = all_folders(~ismember({all_folders.name},{'.','..'}));

sorted_neural_data = [];
uuid_start = 1;

for file_ind_in_dir = 1:length(all_folders)
    filename = all_folders(file_ind_in_dir).name;
    clc;
    disp_str = ['Starting to add file: ' filename ' to sorted neural data| File: ' num2str(file_ind_in_dir) '/' num2str(length(all_folders))];
    disp(disp_str);
    disp('Sorted neural data creation progress:');
    eisg.util.draw_progress_bar(file_ind_in_dir-1, length(all_folders), params.num_ticks_in_progress_bar);
    file_sorting_output = load(fullfile(main_dir, filename, '/matlab/sort.mat'));
    file_spikes = readmda(fullfile(main_dir, filename, '/ms4/firings.mda'));
    file_sorting_output.spikesindices = file_spikes;
    dash_locs = strfind(filename, '-');
    file_sorting_output.region = filename(dash_locs(end)+1:end);
    file_sorting_output = update_sorting_output_struct(file_sorting_output);
    if params.integrate_validity
        unit_validity_in_file = validity{file_ind_in_dir};
        file_sorting_output.validity = unit_validity_in_file;
    end
    n_units = file_sorting_output.n_units;
    uuid = uuid_start:uuid_start+n_units-1;
    uuid_start = uuid_start+n_units;
    file_sorting_output.uuid = uuid;
    file_sorting_output = eisg.preprocess.update_calculated_features_in_file_sorting_output( file_sorting_output, params );
    sorted_neural_data = [sorted_neural_data file_sorting_output];
end
clc;
disp('Sorted neural data creation progress:');
eisg.util.draw_progress_bar(file_ind_in_dir, length(all_folders), params.num_ticks_in_progress_bar);

end

function sorted_data = update_sorting_output_struct(old_sorted_data)

% Reorganize the consolidated data struct

sorted_data                                 = struct();
sorted_data.filename                        = old_sorted_data.src_filename;
sorted_data.region                          = old_sorted_data.region;
sorted_data.n_units                         = length(old_sorted_data.maxchn);
sorted_data.uuid                            = [];
sorted_data.spikeindices                    = old_sorted_data.spikesindices;
% Update maxchn indices for matab
sorted_data.spikeindices(1,:)               = sorted_data.spikeindices(1,:) + 1;
sorted_data.maxchn                          = cell2mat(old_sorted_data.maxchn) + 1;
sorted_data.validity                        = [];
% Have to write a function to prune example waveforms %
sorted_data.example_wfs                     = old_sorted_data.example_wf;
sorted_data.templates                       = prune_ms4_templates(old_sorted_data.templates, old_sorted_data.maxchn);
sorted_data.unit_metrics                    = old_sorted_data.metrics;
sorted_data.template_features_ms4           = old_sorted_data.features;
sorted_data.normalized_templates            = [];
sorted_data.template_timepoints             = [];
sorted_data.calculated_template_features    = [];
sorted_data.preprocess_params               = old_sorted_data.preprocess_params;
sorted_data.sort_params                     = old_sorted_data.sort_params;
sorted_data.postprocess_params              = old_sorted_data.postprocess_params;

% Legend of all the fields in the struct
sorted_data.legend                          = cell(1);
sorted_data.legend{end, 1}                  = 'filename                     -> name of raw file';
sorted_data.legend{end+1, 1}                = 'region                       -> region recorded from';
sorted_data.legend{end+1, 1}                = 'n_units                      -> number of ms4-extracted units';
sorted_data.legend{end+1, 1}                = 'uuid                         -> unique unit id for each sorted unit';
sorted_data.legend{end+1, 1}                = 'spiketindices                -> [maxchn, spikeindex, unitnum] x n_points';
sorted_data.legend{end+1, 1}                = 'maxchn                       -> channel units were extracted from';
sorted_data.legend{end+1, 1}                = 'validity                     -> unit validity after manual curation';
sorted_data.legend{end+1, 1}                = 'example_wfs                  -> 300 or less example spikes per unit';
sorted_data.legend{end+1, 1}                = 'templates                    -> template wfs extracted by ms4/spikeinterface';
sorted_data.legend{end+1, 1}                = 'unit_metrics                 -> firing-related metrics of units';
sorted_data.legend{end+1, 1}                = 'template_features_ms4        -> template shape features retreived from ms4/spikeinterface';
sorted_data.legend{end+1, 1}                = 'normalized_templates         -> norm templates used for feature calculation';
sorted_data.legend{end+1, 1}                = 'template_timepoints          -> timepoints for templates used for feature calculation';
sorted_data.legend{end+1, 1}                = 'calculated_template_features -> template shape features calculated based on extraction from raw';
sorted_data.legend{end+1, 1}                = 'preprocess_params            -> parameters used for preprocessing of raw data';
sorted_data.legend{end+1, 1}                = 'sort_params                  -> parameters used by the ms4 sorter';
sorted_data.legend{end+1, 1}                = 'postprocess_params           -> parameters used for postprocessing of templates';

end


function pruned_template = prune_ms4_templates(templates, maxchn)

maxchn = cell2mat(maxchn) + 1;
pruned_template = nan( numel(maxchn), size(templates, 3) );
for unit = 1:numel(maxchn)
    pruned_template(unit, :) = squeeze( templates(unit, maxchn(unit), :) );
end

end