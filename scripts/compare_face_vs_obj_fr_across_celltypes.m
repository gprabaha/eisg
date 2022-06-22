
clc;
clear;

data_p = fullfile( eisg.util.project_path, 'processed_data');

sorted = shared_utils.io.fload( fullfile(data_p, 'sorted_neural_data_social_gaze.mat') );
events = shared_utils.io.fload( fullfile(data_p, 'events.mat') );

%%
ct_labels = load_cell_type_labels( data_p );

[unit_spike_ts, unit_wfs, spike_labels] = linearize_sorted( sorted );
bfw.add_monk_labels( spike_labels );

[uuid_I, uuids] = findall( spike_labels, 'uuid', find(spike_labels, {'valid-unit', 'maybe-valid-unit'}) );
match_I = bfw.find_combinations( ct_labels, uuids );

for i = 1:numel(uuid_I)
  if ( ~isempty(match_I{i}) )
    ct_label = cellstr( ct_labels, 'cell-type', match_I{i} );
    addsetcat( spike_labels, 'cell-type', ct_label, uuid_I{i} );
  end
end

events = add_whole_face_whole_object_rois( events );
evts = bfw.event_column( events, 'start_time' );

%%
rois = { 'face', 'eyes_nf', 'whole_face', 'right_nonsocial_object_whole_face_matched' };

evt_mask = find( events.labels, [{'m1'}, rois] );
spk_mask = find( spike_labels, {'valid-unit', 'maybe-valid-unit'} );

min_t = -0.5;
max_t = 0.5;
bin_width = 0.01;

[psth_matrix, psth_labels, t] = compute_psth(...
    unit_spike_ts, spike_labels, spk_mask ...
  , evts, events.labels, evt_mask ...
  , min_t, max_t, bin_width );

%%

pre_time  = t>-0.5 & t<=0;
post_time = t>0 & t<=0.5;

spike_pre = nanmean( psth_matrix(:, pre_time), 2 );
spike_post = nanmean( psth_matrix(:, post_time), 2 );

pre_rs_outs = dsp3.ranksum( spike_pre, psth_labels, {'uuid'}, 'whole_face', 'right_nonsocial_object_whole_face_matched' );
post_rs_outs = dsp3.ranksum( spike_post, psth_labels, {'uuid'}, 'whole_face', 'right_nonsocial_object_whole_face_matched' );

% pre_rs_outs = dsp3.ranksum( spike_pre, psth_labels, {'uuid'}, 'face', 'eyes_nf' );
% post_rs_outs = dsp3.ranksum( spike_post, psth_labels, {'uuid'}, 'face', 'eyes_nf' );

% pre_anova_outs = dsp3.anovan( spike_pre, psth_labels, {'uuid'}, {'whole_face', 'right_nonsocial_object_whole_face_matched'} );
% post_anova_outs = dsp3.anovan( spike_post, psth_labels, {'uuid'}, {'whole_face', 'right_nonsocial_object_whole_face_matched'} );

pre_rs_ps = cellfun( @(x) x.p, pre_rs_outs.rs_tables );
post_rs_ps = cellfun( @(x) x.p, post_rs_outs.rs_tables );

pre_rs_zs = cellfun( @(x) x.zval, pre_rs_outs.rs_tables );
post_rs_zs = cellfun( @(x) x.zval, post_rs_outs.rs_tables );

%%

mask = find( psth_labels, {'eyes_nf', 'face', 'right_nonsocial_object_whole_face_matched'} );
[anova_labels, anova_I] = keepeach( psth_labels', 'uuid', mask );
anova_outs = hanova( spike_pre, psth_labels', mask, anova_I );
anova_outs_post = hanova( spike_post, psth_labels', mask, anova_I );
anova_ps = cate1( anova_outs.ps );
anova_post_ps = cate1( anova_outs_post.ps );

%%

factor_ind = 2;% 2: social vs nonsocial; 1: eye vs non-eye face
sig_cells = anova_ps(:, factor_ind) < 0.05; 
sig_cells = sig_cells | anova_post_ps(:, factor_ind) < 0.05;
% [cell_type_I, cell_types] = findall( pre_rs_outs.rs_labels, {'cell-type', 'region'}, findnone(rs_labels, 'b') );
[cell_type_I, cell_types] = findall( anova_labels, {'cell-type', 'region'}, findnone(anova_labels, 'b') );
cell_type_props = cellfun( @(x) pnz(sig_cells(x)), cell_type_I );

prop_table = table( cell_type_props, 'RowNames', fcat.strjoin(cell_types, '_') );
prop_table

%%

chi2_mask = findnone( anova_labels, 'b' );
[count_labels, cell_type_I, cell_types] = keepeach( anova_labels', {'cell-type', 'region'}, chi2_mask );
cell_type_freqs = cellfun( @(x) nnz(sig_cells(x))/numel(x), cell_type_I );

[chi2_info, chi2_labels] = dsp3.chi2_tabular_frequencies( ...
  cell_type_freqs, count_labels, {}, 'cell-type', 'region' );

%%

factor_ind = 2;% 2: social vs nonsocial; 1: eye vs non-eye face

use_labels = anova_labels';
sig_cells = anova_ps(:, factor_ind) < 0.05; 
sig_cells = sig_cells | anova_post_ps(:, factor_ind) < 0.05;

chi2_mask = findnone( anova_labels, 'b' );
[count_labels, cell_type_I, cell_types] = keepeach( anova_labels', {'cell-type', 'region'}, chi2_mask );
cell_type_freqs = cellfun( @(x) nnz(sig_cells(x))/numel(x), cell_type_I );

% cell_type_props = cellfun( @(x) sum(sig_cells(x)) / numel(sig_cells(x)), cell_type_I );
cell_type_props = cellfun( @(x) numel(sig_cells(x)), cell_type_I );
prop_table = table( cell_type_props, 'RowNames', fcat.strjoin(cell_types, '_') );
prop_table

[tbl, reg_labels] = do_prop_tests( sig_cells, use_labels', findnone(use_labels, 'b') );
tbl

cell_type_props = cellfun( @(x) sum(sig_cells(x) & post_rs_zs(x) > 0) / sum(sig_cells(x)), cell_type_I );
% cell_type_props = cellfun( @(x) sum(sig_cells(x) & pre_rs_zs(x) > 0) / sum(sig_cells(x)), cell_type_I );
prop_table = table( cell_type_props, 'RowNames', fcat.strjoin(cell_types, '_') );
prop_table

post_face_high_cells = post_rs_zs > 0;
% post_face_high_cells = pre_rs_zs > 0;
use_labels = fcat.from( anova_labels( find(sig_cells), : ), getcats( anova_labels ) )';
[tbl, reg_labels] = do_prop_tests( post_face_high_cells( sig_cells ), use_labels', findnone(use_labels, 'b') );
tbl

%%

factor_ind = 1;
sig_eye_cells = sig_cells & ( anova_ps(:, factor_ind) < 0.05 | anova_post_ps(:, factor_ind) < 0.05 );
sig_eye_cells = sig_eye_cells(sig_cells);
use_labels = fcat.from( anova_labels( find(sig_cells), : ), getcats( anova_labels ) )';

[count_labels, cell_type_I, cell_types] = keepeach( use_labels', {'cell-type', 'region'}, findnone(use_labels, 'b') );
%cell_type_props = cellfun( @(x) sum(sig_eye_cells(x)) / numel(sig_eye_cells(x)), cell_type_I );
cell_type_props = cellfun( @(x) sum(sig_eye_cells(x)), cell_type_I );
prop_table = table( cell_type_props, 'RowNames', fcat.strjoin(cell_types, '_') );
prop_table

[tbl, reg_labels] = do_prop_tests( sig_eye_cells, use_labels', findnone(use_labels, 'b') );
tbl
%%

rs_labels = pre_rs_outs.rs_labels';

% sig_cells = pre_rs_ps < 0.05 | post_rs_ps < 0.05;
% [cell_type_I, cell_types] = findall( pre_rs_outs.rs_labels, {'cell-type', 'region'}, findnone(rs_labels, 'b') );
[cell_type_I, cell_types] = findall( pre_rs_outs.rs_labels, {'cell-type', 'region'}, findnone(rs_labels, 'b') );
cell_type_props = cellfun( @(x) sum(sig_cells(x) & post_rs_zs(x) > 0) / sum(sig_cells(x)), cell_type_I );

prop_table = table( cell_type_props, 'RowNames', fcat.strjoin(cell_types, '_') );
prop_table

%%

function [tbl, reg_labels, tbls] = do_prop_tests(meets_criterion, labels, mask)

assert_ispair( meets_criterion, labels );

[reg_labels, reg_I] = keepeach( labels', 'region', mask );
n_I = cellfun( @(x) find(labels, 'n', x), reg_I, 'un', 0 );
m_I = cellfun( @(x) find(labels, 'm', x), reg_I, 'un', 0 );

tbls = do_prop_test( meets_criterion, n_I, m_I );
chi2ps = [ tbls.p ]';
tbl = array2table( chi2ps, 'RowNames', fcat.strjoin(reg_labels(:, 'region')') );

end

function tbls = do_prop_test(meets_criterion, inds_a, inds_b)

xa = cellfun( @(x) sum(meets_criterion(x)), inds_a );
xb = cellfun( @(x) sum(meets_criterion(x)), inds_b );
na = cellfun( @numel, inds_a );
nb = cellfun( @numel, inds_b );

tbls = struct( 'h', {}, 'p', {}, 'chi2stat', {} );
for i = 1:numel(xa)
  [h, p, chi2stat] = prop_test( [xa(i), xb(i)], [na(i), nb(i)], false );
  tbls(i) = struct( 'h', h, 'p', p, 'chi2stat', chi2stat );
end

end


