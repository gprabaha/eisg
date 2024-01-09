function plot_zscored_auc_timeseries_across_celltypes(...
    t, z_scored_aucs, auc_labels, region, excluded_categories)

unit_mask = pipe( rowmask( auc_labels ), ...
    @(m) find( auc_labels, region, m ), ...
    @(m) findnone( auc_labels, excluded_categories, m )...
    );

abs_zscored_auc = abs(z_scored_aucs);

custom_narrow_map = [linspace(1, 0, 100)', linspace(1, 0.3294, 100)', linspace(1, 0.5216, 100)']; % Kinda blue #005485
custom_broad_map = [linspace(1, 0.7804, 100)', linspace(1, 0.3216, 100)', linspace(1, 0.1647, 100)']; % Kinda orange #c7522a

[I, C] = findall( auc_labels, {'cell-type'}, unit_mask );
figure();
axs = plots.panels( numel(I) );

for i = 1:numel(axs)
  ind = I{i};
  abs_zscore_auc_subset = abs_zscored_auc(ind,:);
  [~, max_auc_loc] = max( abs_zscore_auc_subset, [], 2);
  [~,sorted_ind] = sort( max_auc_loc );
  imagesc( axs(i), t, 1:numel(ind), abs_zscore_auc_subset(sorted_ind, :) );
  colorbar( axs(i) );
  set( axs(i), 'clim', [0, max(max(abs(z_scored_aucs)))] );
  if strcmp(C{i}, 'narrow')
      custom_map = custom_broad_map;
  elseif strcmp(C{i}, 'broad')
      custom_map = custom_narrow_map;
  else
      error_message = sprintf('No colormap for celltype: %s', char(C{i}));
      error(error_message);
  end
  colormap( axs(i), custom_map );
  title( axs(i), strrep(fcat.strjoin(C(:, i), ' | '), '_', ' ') );
end
set(gcf, 'Position',  [200, 200, 900, 500]);

end 

