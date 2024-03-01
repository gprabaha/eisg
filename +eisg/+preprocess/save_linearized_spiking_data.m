function save_linearized_spiking_data()

data_p = fullfile( eisg.util.project_path, 'processed_data' );

disp('Loading data...');
sorted = shared_utils.io.fload( fullfile(data_p,...
  'sorted_neural_data_social_gaze.mat') );
ct_labels = shared_utils.io.fload(fullfile(data_p,...
    'celltype-labels_pfc-combined-class_p2v.mat'), 'ct_labels');
disp('Done');

disp('Processing data...');
validity_filter = {'valid-unit', 'maybe-valid-unit'};
[all_unit_spike_ts, all_unit_wfs, unit_labels] = eisg.util.linearize_sorted(sorted);
bfw.add_monk_labels(unit_labels);
[uuid_I, uuids] = findall(unit_labels, 'uuid', find(unit_labels, validity_filter));
match_I = bfw.find_combinations(ct_labels, uuids);
for i = 1:numel(uuid_I)
    if (~isempty(match_I{i}))
        ct_label = cellstr(ct_labels, 'cell-type', match_I{i});
        addsetcat(unit_labels, 'cell-type', ct_label, uuid_I{i});
    end
end
replace(unit_labels, 'n', 'narrow');
replace(unit_labels, 'm', 'broad');
replace(unit_labels, 'b', 'outlier');
disp('Done');

spike_data_filename     = 'spike_data_celltype_labelled.mat';
savepath = fullfile( data_p, spike_data_filename );
disp(['Saving: ', savepath ]);
save( savepath, ...
    'all_unit_spike_ts',  'all_unit_wfs', 'unit_labels');
disp('Done')

end