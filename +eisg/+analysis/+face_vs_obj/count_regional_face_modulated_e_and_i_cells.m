function [counting_mat, counting_mat_cols] = count_regional_face_modulated_e_and_i_cells(valid_unit_psth, varargin)

defaults = eisg.util.make_analysis_params_struct();
defaults.tasktype = 'free_viewing';
defaults.pre_fix_time_win = [-0.4 -0.05];
defaults.post_fix_time_win = [0.05 0.4];

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
% Add celltype category to psth label data
addcat(all_psth_labels, 'celltype');
regions = combs(all_psth_labels, 'region');
unit_counter = 0;
total_valid_units = 0;
for region = regions
    ei_labels = eisg.util.fetch_ei_labels( region, params );
    label_mat = ei_labels.label_mat;
    label_cols = ei_labels.label_mat_cols;
    uuid_col_in_ei_mat = strcmp( label_cols, 'uuid' );
    ei_label_col_in_ei_mat = strcmp( label_cols, 'ei_label' );
    regional_uuids = label_mat(:, uuid_col_in_ei_mat)';
    n_units = length(regional_uuids);
    total_valid_units = total_valid_units + n_units;
    for uuid_index = 1:n_units
        uuid = label_mat(uuid_index, uuid_col_in_ei_mat);
        uuid_ind_in_psth = find( all_psth_labels, char(uuid) );
        uuid_ei_label = label_mat(uuid_index, ei_label_col_in_ei_mat);
        setcat( all_psth_labels, 'celltype', cellstr(uuid_ei_label) , uuid_ind_in_psth );
    end
end
addcat(all_psth_labels, 'roi_simp');
face_ind = find( all_psth_labels, face_roi );
setcat( all_psth_labels, 'roi_simp', {'in_face'}, face_ind );
obj_ind = find( all_psth_labels, obj_roi );
setcat( all_psth_labels, 'roi_simp', {'in_obj'}, obj_ind );
tasktype_inds = find(all_psth_labels, {tasktype, 'in_face', 'in_obj'});
all_psth = all_psth(tasktype_inds, :);
all_psth_labels = all_psth_labels(tasktype_inds);

prune( all_psth_labels );

counting_mat = [];
for region = regions
    regional_inds = find(all_psth_labels, region);
    [regional_unit_inds, regional_units] = findall( all_psth_labels(regional_inds), 'uuid' );
    for regional_unit_num = 1:numel(regional_units)
        clc;
        disp('Progress:');
        eisg.util.draw_progress_bar(unit_counter, total_valid_units, params.num_ticks_in_progress_bar);
        unit_counter = unit_counter+1;
        unit_inds_in_psth = regional_unit_inds{regional_unit_num};
        face_activity = all_psth(find( all_psth_labels(unit_inds_in_psth), 'in_face' ), :);
        obj_activity = all_psth(find( all_psth_labels(unit_inds_in_psth), 'in_obj' ), :);
        sig_face = is_unit_sig(face_activity, obj_activity, bin_t, params);
        if sig_face
            eye_activity = all_psth(find( all_psth_labels(unit_inds_in_psth), 'eye' ), :);
            ne_face_activity = all_psth(find( all_psth_labels(unit_inds_in_psth), 'ne_face' ), :);
            sig_eye = is_unit_sig(eye_activity, ne_face_activity, bin_t, params);
        else
            sig_eye = 0;
        end
        celltype = all_psth_labels(unit_inds_in_psth(1), 'celltype');
        cell_prop = categorical([region, regional_units(regional_unit_num), celltype, {num2str(sig_face)}, {num2str(sig_eye)}]);
        counting_mat = [counting_mat; cell_prop];
    end
end

counting_mat_cols = {'region', 'uuid', 'celltype', 'sig_face', 'sig_eye'};

end

function tf = is_unit_sig(roi_activity, ctrl_activity, bin_t, params)

% Pre-window analysis
pre_inds = bin_t>params.pre_fix_time_win(1) & bin_t<params.pre_fix_time_win(2);
post_inds = bin_t>params.post_fix_time_win(1) & bin_t<params.post_fix_time_win(2);

roi_pre = remove_nan_rows( roi_activity(:,pre_inds) );
roi_post = remove_nan_rows( roi_activity(:,post_inds) );
ctrl_pre = remove_nan_rows( ctrl_activity(:,pre_inds) );
ctrl_post = remove_nan_rows( ctrl_activity(:,post_inds) );

pre_roi_sig = 0;
if ~( isempty(roi_pre) || isempty(ctrl_pre) )
    pre_roi_sig = ranksum_compare_activities( roi_pre, ctrl_pre );
end
post_roi_sig = 0;
if ~( isempty(roi_post) || isempty(ctrl_post) )
    post_roi_sig = ranksum_compare_activities( roi_post, ctrl_post );
end

tf = pre_roi_sig | post_roi_sig;

end

function h = ranksum_compare_activities( mat_1, mat_2 )

mean_mat_1 = nanmean(mat_1, 2);
mean_mat_2 = nanmean(mat_2, 2);
[~,h] = ranksum(mean_mat_1, mean_mat_2);

end

function new_mat = remove_nan_rows(old_mat)

new_mat = old_mat(~any( ismissing( old_mat ), 2 ), :);

end

