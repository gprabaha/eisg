% Calculate probability estimates and confidence intervals for 'broad' and 'narrow' within 'acc'
[phat_broad_acc, pci_broad_acc] = binofit(sum(broad_acc_data), numel(broad_acc_data));
[phat_narrow_acc, pci_narrow_acc] = binofit(sum(narrow_acc_data), numel(narrow_acc_data));

% Calculate probability estimates and confidence intervals for 'broad' and 'narrow' within 'bla'
[phat_broad_bla, pci_broad_bla] = binofit(sum(broad_bla_data), numel(broad_bla_data));
[phat_narrow_bla, pci_narrow_bla] = binofit(sum(narrow_bla_data), numel(narrow_bla_data));

% Display probability estimates and confidence intervals for 'broad' and 'narrow' within 'acc'
disp('Broad within Acc:');
disp(['Probability Estimate: ', num2str(phat_broad_acc)]);
disp(['Confidence Intervals: [', num2str(pci_broad_acc(1)), ', ', num2str(pci_broad_acc(2)), ']']);

disp('Narrow within Acc:');
disp(['Probability Estimate: ', num2str(phat_narrow_acc)]);
disp(['Confidence Intervals: [', num2str(pci_narrow_acc(1)), ', ', num2str(pci_narrow_acc(2)), ']']);

% Display probability estimates and confidence intervals for 'broad' and 'narrow' within 'bla'
disp('Broad within Bla:');
disp(['Probability Estimate: ', num2str(phat_broad_bla)]);
disp(['Confidence Intervals: [', num2str(pci_broad_bla(1)), ', ', num2str(pci_broad_bla(2)), ']']);

disp('Narrow within Bla:');
disp(['Probability Estimate: ', num2str(phat_narrow_bla)]);
disp(['Confidence Intervals: [', num2str(pci_narrow_bla(1)), ', ', num2str(pci_narrow_bla(2)), ']']);

% Calculate the proportion of overlap for 'broad' and 'narrow' within 'acc'
overlap_acc = calculate_overlap(pci_broad_acc, pci_narrow_acc);

% Calculate the proportion of overlap for 'broad' and 'narrow' within 'bla'
overlap_bla = calculate_overlap(pci_broad_bla, pci_narrow_bla);

% Display the proportion of overlap
disp(['Proportion of Overlap for ''broad'' and ''narrow'' within ''acc'': ', num2str(overlap_acc)]);
disp(['Proportion of Overlap for ''broad'' and ''narrow'' within ''bla'': ', num2str(overlap_bla)]);

%%

% Import the Statistics and Machine Learning Toolbox if not already done
if ~license('test', 'Statistics_Toolbox')
    error('Statistics and Machine Learning Toolbox is required for this analysis.');
end
% Perform Fisher's Exact Test for 'broad' and 'narrow' within 'acc'
[~, p_value_acc] = fishertest([sum(broad_acc_data), numel(broad_acc_data) - sum(broad_acc_data); ...
                                sum(narrow_acc_data), numel(narrow_acc_data) - sum(narrow_acc_data)]);

% Perform Fisher's Exact Test for 'broad' and 'narrow' within 'bla'
[~, p_value_bla] = fishertest([sum(broad_bla_data), numel(broad_bla_data) - sum(broad_bla_data); ...
                                sum(narrow_bla_data), numel(narrow_bla_data) - sum(narrow_bla_data)]);

% Display the p-values
disp(['Fisher''s Exact Test p-value for ''broad'' and ''narrow'' within ''acc'': ', num2str(p_value_acc)]);
disp(['Fisher''s Exact Test p-value for ''broad'' and ''narrow'' within ''bla'': ', num2str(p_value_bla)]);


%%

% Helper function to calculate overlap proportion
function overlap_proportion = calculate_overlap(ci1, ci2)
    overlap_start = max(ci1(1), ci2(1));
    overlap_end = min(ci1(2), ci2(2));
    
    if overlap_end < overlap_start
        overlap_proportion = 0;
    else
        total_span = max(ci1(2), ci2(2)) - min(ci1(1), ci2(1));
        overlap_proportion = (overlap_end - overlap_start) / total_span;
    end
end
