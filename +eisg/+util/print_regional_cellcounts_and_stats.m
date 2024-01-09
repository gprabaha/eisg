function print_regional_cellcounts_and_stats(counting_mat, counting_mat_cols)

region_col = strcmp(counting_mat_cols, 'region');
uuid_col = strcmp(counting_mat_cols, 'uuid');
celltype_col = strcmp(counting_mat_cols, 'celltype');
face_col = strcmp(counting_mat_cols, 'sig_face');
eye_col = strcmp(counting_mat_cols, 'sig_eye');


regions = unique( counting_mat(:, region_col) );
summary = [];
for region = regions'
    region_inds = counting_mat(:, region_col) == region;
    regional_mat = counting_mat(region_inds, :);
    unit_count = sum(region_inds);
    celltypes = unique( regional_mat(:, celltype_col) );
    for celltype = celltypes'
        celltype_inds = regional_mat(:, celltype_col) == celltype;
        celltype_mat = regional_mat(celltype_inds, :);
        celltype_count = sum(celltype_inds);
        face_sig_unit_inds = celltype_mat(:, face_col) == '1';
        face_sig_unit_count = sum(face_sig_unit_inds);
        eye_sig_unit_inds = celltype_mat(:, eye_col) == '1';
        eye_sig_unit_count = sum(eye_sig_unit_inds);
        summary_cell = {char(region), unit_count, char(celltype), celltype_count, face_sig_unit_count, eye_sig_unit_count};
        summary = [summary; summary_cell];
    end
    fprintf('In region %s there are %d units\n', summary{end,1}, summary{end,2} );
    n_e_cells = summary{end-1,4};
    n_e_face_cells = summary{end-1,5};
    n_i_cells = summary{end,4};
    n_i_face_cells = summary{end,5};
    p = eisg.util.prop_test([n_e_face_cells, n_i_face_cells], [n_e_cells n_i_cells], false);
    fprintf('Face activated E cells/Total E cells: %d/%d; Face activated I cells/Total I cells: %d/%d; Chi-sq p: %0.3f\n', n_e_face_cells, n_e_cells, n_i_face_cells, n_i_cells, p);
    n_face_cells = n_e_face_cells+n_i_face_cells;
    p = eisg.util.prop_test([n_e_face_cells, n_i_face_cells], [n_face_cells n_face_cells], false);
    fprintf('Face activated E cells/Total face cells: %d/%d; Face activated I cells/Total face cells: %d/%d; Chi-sq p: %0.3f\n', n_e_cells, n_e_face_cells, n_face_cells, n_face_cells, p);
    
end


end