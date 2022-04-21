function [aic,bic,fom_aic,fom_bic,nclus_aic,nclus_bic] = AICBICAnalysis(datanonnorm4Cluster, labels, clustFilt, minNumOfCluster, maxNumOfCluster, offset)
% This function computes the Akaike information criterion and Bayesian
% information criterion for every k input to the kMeansClustering. The
% elbow point of AIC and BIC curve is the optimal number of clusters.
% Author: Ethan Wang
% Date: 09-01-2021
% Inputs:
% datanonnorm4Cluster: non-normalized data for clustering
% labels: cluster id for each unit
% clustFilt: unit id for every cluster
% minNumOfCluster
% maxNumOfCluster
% offset = min(k) - 1
% Output:
% aic
% bic
% fom_aic: the ratios of aic changes for adjacent k
% fom_bic: the ratios of bic changes for adjacent k
% nclus_aic: the elbow of the aic curve
% nclus_bic: the elbow of the bic curve


minNumOfElementsCluster = 5;
for k = minNumOfCluster:maxNumOfCluster % num of clusters
    relIndex = k-offset;
    %control.data = tmp.datanorm4Cluster;
    % zscore
    control.data = zscore(datanonnorm4Cluster);
    numDim = size(control.data,2);
    nElements = length(labels);
    control.labels = nan(1,nElements);
    control.clustOrder = 1:max(labels(relIndex,:));
    cont = 0; % count of real clusters for each k
    for i = control.clustOrder % exlcuding clusters with less than 5 units
        if (length(clustFilt{relIndex,i}) < minNumOfElementsCluster)
         continue
        else
          cont = cont + 1; 
          control.labels(labels(relIndex,:)==i) = cont; % record its cluster id
          control.centers(cont,:) = mean(control.data(control.labels==cont,:),1); % record the cluster's mean value
          control.elemInCluster{cont} = clustFilt{relIndex,i}; % record the cluster's units ids
          wcd{cont} = control.data(control.labels==cont,:) - control.centers(cont,:); % within cluster difference
          withinss(cont) = sumsqr(wcd{cont}); % within cluster sum of squares
        end
    end
    for i = 1:cont
        % number of units in each real cluster
        Nc(i) = length(control.elemInCluster{i}); 
        % variance of each feature within each cluster normalizing by n-1
        Vc(:,i) = var(control.data(control.labels==i,:),0,1);
    end
    %clear tmp control;

    numRelClusters = cont;
    
    %logL = -Nc.*sum(.5*log(Vc),1); % there is no cluster with only 1 element
    %bic(relIndex) = sum(Vc)+numRelClusters*numDim*log(nElements);
    
    % method 1
    % aic(relIndex) = withinss + 2*numDim*numRelClusters;
    % bic(relIndex) = withinss + log(nElements)*numDim*numRelClusters;
    
    % method 2 (requires zscored data)
    % variance of each feature across all units
    V = var(control.data,0,1)'*ones(1,numRelClusters);
    %logL = -Nc.*sum(.5*log(Vc+V),1);
    logL = -Nc.*log(withinss+1e-4);
    aic(relIndex) = -2 * sum(logL) + 4*numRelClusters*numDim;
    bic(relIndex) = -2 * sum(logL) + 2*numRelClusters*numDim*log(nElements);
end

% find the elbow
aic_norm = (aic-min(aic))/(max(aic)-min(aic));
bic_norm = (bic-min(bic))/(max(bic)-min(bic));
% v_aic = -diff(aic);
% v_bic = -diff(bic);
v_aic = -diff(aic_norm);
v_bic = -diff(bic_norm);
nv = length(v_aic);
fom_aic = v_aic(1:(nv-1))./v_aic(2:nv);
fom_bic = v_bic(1:(nv-1))./v_bic(2:nv);
% the elbow is defined as the point where the ratio of changes is the biggest
nclus_aic = find(fom_aic==max(fom_aic),1) + 1; 
nclus_bic = find(fom_bic==max(fom_bic),1) + 1;