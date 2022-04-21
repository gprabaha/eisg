function plot_psth_for_each_unit( valid_unit_psth, varargin )

defaults = eisg.util.make_analysis_params_struct();
defaults.tasktype = 'free_viewing';
defaults.face_roi = {'face', 'eyes_nf'};
defaults.eye_roi = {'eyes_nf'};
defaults.ne_face_roi = {'face'};
defaults.obj_roi = {'left_nonsocial_object', 'right_nonsocial_object'};

params = shared_utils.general.parsestruct( defaults, varargin );

all_psth = valid_unit_psth.trial_psth;
bin_t = valid_unit_psth.bin_t;
bin_width = valid_unit_psth.bin_width;
[~, uuid_ind] = ismember( {'uuid'}, valid_unit_psth.label_categories );
all_psth_labels = fcat.from( ...
     valid_unit_psth.trial_psth_labels(:, 1:uuid_ind) ...
    , valid_unit_psth.label_categories(1:uuid_ind) );


face_roi = params.face_roi;
eye_roi = params.eye_roi;
ne_face_roi = params.ne_face_roi;
obj_roi = params.obj_roi;
tasktype = params.tasktype;

% Pruning out the labels
addcat(all_psth_labels, 'roi_simp');
%
face_ind = find( all_psth_labels, face_roi );
setcat( all_psth_labels, 'roi_simp', {'in_face'}, face_ind );
obj_ind = find( all_psth_labels, obj_roi );
setcat( all_psth_labels, 'roi_simp', {'in_obj'}, obj_ind );

% eye_ind = find( all_psth_labels, eye_roi );
% setcat( all_psth_labels, 'roi_simp', {'in_eyes'}, eye_ind );
% ne_face_ind = find( all_psth_labels, ne_face_roi );
% setcat( all_psth_labels, 'roi_simp', {'in_ne_face'}, ne_face_ind );

% Remove all epochs other than rwdOn and inTarget
tasktype_inds = find(all_psth_labels, {tasktype, 'in_face', 'in_obj'});
all_psth = all_psth(tasktype_inds, :);
all_psth_labels = all_psth_labels(tasktype_inds);
prune( all_psth_labels );

regions = combs(all_psth_labels, 'region');

unit_counter = 0;
total_valid_units = 0;
for region = regions
    ei_labels = eisg.util.fetch_ei_labels( region, params );
    label_mat = ei_labels.label_mat;
    n_units = size(label_mat, 1);
    total_valid_units = total_valid_units + n_units;
end
for region = regions
    ei_labels = eisg.util.fetch_ei_labels( region, params );
    label_mat = ei_labels.label_mat;
    label_cols = ei_labels.label_mat_cols;
    uuid_col_in_ei_mat = strcmp( label_cols, 'uuid' );
    ei_label_col_in_ei_mat = strcmp( label_cols, 'ei_label' );
    regional_uuids = label_mat(:, uuid_col_in_ei_mat)';
    for uuid_index = 1:length(regional_uuids)
        clc;
        disp(['Unit PSTH plotting progress (in ' char(region) '):']);
        eisg.util.draw_progress_bar(unit_counter, total_valid_units, params.num_ticks_in_progress_bar);
        unit_counter = unit_counter + 1;
        uuid = label_mat(uuid_index, uuid_col_in_ei_mat);
        uuid_ei_label = label_mat(uuid_index, ei_label_col_in_ei_mat);
        uuid_inds_in_psth = find(all_psth_labels, char(uuid));
        uuid_psth = all_psth(uuid_inds_in_psth, :);
        uuid_psth_labels = all_psth_labels(uuid_inds_in_psth);
        pl = plotlabeled.make_common();
        pl.per_panel_labels = true;
        pl.one_legend = false;
        pl.x = bin_t;
        pl.summary_func = @nanmean;
        pl.add_errors = false;
        pl.error_func = @plotlabeled.nansem;
        pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.7);
        pl.add_smoothing = true;
        pl.main_line_width = 1.5;
        pl.error_line_width = 0.5;
        pl.lines( uuid_psth./bin_width, uuid_psth_labels, {'roi_simp'}, {'region', 'uuid'} );
        
        % Saving fig
        save_folder_path = fullfile( './', params.plot_folder, 'single_unit_psth', char(region), char(uuid_ei_label) );
        if ~exist(save_folder_path, 'dir')
            mkdir(save_folder_path)
        end
        dsp3.req_savefig( gcf, save_folder_path, uuid_psth_labels, {'region', 'uuid'}, 'single_unit_psth' );
        
    end
end
clc;
disp('Unit PSTH plotting progress:');
eisg.util.draw_progress_bar(unit_counter, total_valid_units, params.num_ticks_in_progress_bar);


end