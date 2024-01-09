function [unit_wf_features, feature_labels] = append_mean_fr_in_wf_features( unit_wf_features, feature_labels, unit_mean_fr, mean_fr_labs )

unit_wf_features = [unit_wf_features; unit_mean_fr(:,1)];
addsetcat( mean_fr_labs, 'wf_feature', 'mean_fr' );
feature_label_cats = getcats(feature_labels);
fr_label_subset = fcat.from( mean_fr_labs(:, feature_label_cats'), feature_label_cats );
append( feature_labels, fr_label_subset );

end