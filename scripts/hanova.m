function out = hanova(data, labels, mask, anova_I)

assert_ispair( data, labels );

social_group = make_social_group( labels, mask );
specified_factors = { 'roi', 'social' };
anova_factors = { 'roi' };

ps = cell( size(anova_I) );
stat_tables = cell( size(ps) );

for i = 1:numel(anova_I)  
  anova_ind = anova_I{i}; 
  
  subset_spikes = data(anova_ind);
  subset_social_group = social_group(anova_ind);
  
  groups = cellfun( @(x) categorical(labels, x, anova_ind), anova_factors, 'un', 0 );
  
  if ( ismember('social', specified_factors) )
    groups{end+1} = subset_social_group;
  end
  
  nesting = zeros( numel(anova_factors)+1 );
  roi_ind = find( strcmp(anova_factors, 'roi') );
  % Nest rois in social
  nesting(roi_ind, numel(groups)) = 1;

  [p, tbl, stats] = anovan( subset_spikes, groups ...
    , 'nested', nesting ...
    , 'varnames', [anova_factors, {'social'}] ...
    , 'display', 'off' ...
  );

  ps{i} = p(:)';
  stat_tables{i} = stats;
  
%   social_rois = sort( combs(labels, 'roi', anova_ind(subset_social_group == 'social')) );
%   nonsocial_rois = sort( combs(labels, 'roi', anova_ind(subset_social_group == 'nonsocial')) );
%   [roi_I, roi_combs] = findall( labels, 'roi', anova_ind );
%   
%   [all_rois, sort_ind] = sort( roi_combs );
%   roi_I = roi_I(sort_ind);
%   
%   if ( ismember('social', specified_factors) )
%     post_hoc_groups = { social_rois, nonsocial_rois, all_rois };
%     post_hoc_group_names = { 'social', 'nonsocial', 'all' };
%   else
%     post_hoc_groups = { all_rois };
%     post_hoc_group_names = { 'all' };
%   end
%   
%   post_hoc_info = struct();
%   
%   for j = 1:numel(post_hoc_groups)
%     [post_hoc_p, roi_strs] = post_hoc_comparisons( spikes, labels, post_hoc_groups{j}, anova_ind, post_hoc_func );
%     post_hoc_info.(post_hoc_group_names{j}).p = post_hoc_p;
%     post_hoc_info.(post_hoc_group_names{j}).rois = roi_strs;
%   end
%   
%   group_means = bfw.row_nanmean( spikes, roi_I );
%   [~, group_ordering] = sort( group_means );
%   
%   group_ordering_ids{i} = strjoin( arrayfun(@num2str, group_ordering, 'un', 0), ',' );
%   all_roi_combs{i} = roi_combs;
%   
%   all_post_hoc_info{i} = post_hoc_info;
end

out = struct();
out.stat_tables = stat_tables;
out.ps = ps;

end

function social_group = make_social_group(labels, mask)

social_ind = findor( labels, social_rois(), mask );  
nonsocial_ind = findor( labels, nonsocial_rois(), mask );
social_group = categorical();
social_group(rowmask(labels), 1) = '<undefined>';
social_group(social_ind) = 'social';
social_group(nonsocial_ind) = 'nonsocial';

end

function rois = nonsocial_rois()

rois = { 'left_nonsocial_object', 'right_nonsocial_object', 'nonsocial_object' ...
  , 'left_nonsocial_object_eyes_nf_matched', 'right_nonsocial_object_eyes_nf_matched' ...
  , 'nonsocial_object_eyes_nf_matched' ...
  , 'right_nonsocial_object_whole_face_matched' ...
  , 'left_nonsocial_object_whole_face_matched' ...
};

end

function rois = social_rois()

rois = { 'eyes_nf', 'face', 'face_non_eyes' };

end