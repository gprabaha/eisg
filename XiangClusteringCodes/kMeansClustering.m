function [clustFilt,labels,percElements,centroids] = kMeansClustering(data4Cluster,minNumOfClutser,maxNumOfCluster)


% This function runs different realizations of the Kmeans/Kmedians clustering
% alogorithm with a previously established number of clusters to evaluate
% to what extend clusters are composed of the same elements
% Author: Ethan Wang
% Date: 2021-06-03
% Input:
% datanorm4Cluster: normalized data for cluster, neuron x feature
% minNumOfClutser
% maxNumOfClutser
% Output:
% clustFilt: unit ID within each cluster
% labels: label of clusters for each unit
% percElements: percentage of classified units
% centroids: centroids for every cluster

% to control the time of the analysis
timeId = tic;

% normalize data
data4Cluster = standardization(data4Cluster);

% parameters
numberOfRepeats = 500;   % number of realizations
numOfReplicates = 50;    % number of replicates

% number of clusters to consider (based on
% main_estimate_numberOfClusters.m)
numC = minNumOfClutser:maxNumOfCluster;
cutoff_clusterElements = 5;   % at least 4 or 5 or 6 elements in a cluster to be considered

% number of elements (neurons, in rows) and dimensions (measures, in
% columns)
[nrow,~] = size(data4Cluster);

% initializing labels (clusters to which neurons belong to)
labels = zeros(length(numC),nrow);

ind = 0;
for numClusters = numC
    ind = ind +1;
    fprintf('\n>>>> starting clustering algorithm with %d clusters', numClusters);
    
    pairwiseCluster = tril(nan(nrow,nrow),0) + triu(zeros(nrow),1); % initializing the matrix of pariwise same cluster belonging
    classlabel = zeros(nrow,numberOfRepeats); % it keeps record of neurons' cluster through realizations
    
    numOfIntersections = zeros(numClusters,numClusters);
    minNumOfIntersections = nrow*ones(numClusters,1);
    maxNumOfElem = zeros(numClusters,1);
    
    for i = 1:numberOfRepeats
        
        % part 1: running clustering algorithm
        
        if (mod(i,50)==0)
            fprintf('\n++++ running clustering alogrithm, iteration=%d',i);
        end
        
        [classlabel(:,i),centroidx] = kmeans(data4Cluster,numClusters,'replicates',numOfReplicates,'emptyaction','singleton','distance','sqEuclidean');
        
        if(i>1)
            Qpre = Q;
        end
        Q = ind2cluster(classlabel(:,i));
        
        % part 2: control for neurons belonging to same clusters
        
        ns = [];
        for j = 1:numel(Q)
            ns(j) = numel(Q{j});
            for k = 1:ns(j)
                for l = k+1:ns(j)
                    pairwiseCluster(Q{j}(k),Q{j}(l)) = nansum(pairwiseCluster(Q{j}(k),Q{j}(l))+1);
                end
            end
        end
    end
    
    pairwiseCluster = pairwiseCluster/numberOfRepeats;
    symPairwiseCluster = tril(pairwiseCluster')+triu(pairwiseCluster);
    
    cutoff = 0.9;
    symPairwiseClusterMod = zeros(size(pairwiseCluster));
    symPairwiseClusterMod(symPairwiseCluster>=cutoff) = symPairwiseCluster(symPairwiseCluster>=cutoff);
    reorderedSymPairwiseCluster = symrcm(symPairwiseClusterMod);
    
    pairwiseClusterReordered = triu(symPairwiseClusterMod(reorderedSymPairwiseCluster,reorderedSymPairwiseCluster));
    
    k = 0;
    allElements = [];
    for i = 1:nrow
        rowElements = find(pairwiseClusterReordered(i,:)>=cutoff);
        if(~isempty(rowElements))
            rowElements = union(rowElements,i);
            if(isempty(intersect(rowElements,allElements)))
                k = k+1;
                clust{k} = rowElements;
            else
                clust{k} = union(clust{k},rowElements);
            end
            allElements = union(allElements,rowElements);
        end
    end
    
    % check again that all clusters don't share elements
    flag = 1;
    while flag == 1
        flag = 0;
        l = length(clust);
        for i = 1:(length(clust)-1)
            if(isempty(clust{i}))
                clust{i} = clust{i+1};
                clust{i+1} = {};
            elseif(~isempty(intersect(clust{i},clust{i+1})))
                clust{i} = union(clust{i},clust{i+1});
                clust{i+1} = {};
                l = l-1;
                flag = 1;
            end
        end
        clust = clust(1:l);
    end
    
    cont = 0;
    allElements = [];
    for i = 1:length(clust)
        if(numel(clust{i})>=cutoff_clusterElements)
            cont = cont + 1;
            for j = 1:numel(clust{i})
                clustFilt{ind,cont}(j) = reorderedSymPairwiseCluster(clust{i}(j));
                allElements = [allElements,clustFilt{ind,cont}(j)];
            end
            clustFilt{ind,cont} = sort(clustFilt{ind,cont});
        end
    end
    clear clust;
    
    percElem = 100*length(allElements)/nrow;
    format compact
    display(percElem); 
    percElements(ind) = 100*length(allElements)/nrow;
    
    labels(ind,:) = cluster2label(clustFilt(ind,:),nrow);
end

% compute the centroids for each cluster for each k
for n = 1:size(labels,1)
    for c = 1:max(labels(n,:))
        clear dataInCluster
        dataInCluster = data4Cluster(labels(n,:)==c,:);
        centroids{n,c} = mean(dataInCluster);
    end
end

% total elapsed time in the analysis
elapsedTime = toc(timeId);

elapsedTimeHours = floor(elapsedTime/3600);
elapsedTimeMinuts = floor(mod(elapsedTime,3600)/60);
elapsedTimeSeconds = mod(mod(elapsedTime,3600),60);

fprintf('\n>>>> total time = %d hours, %d minutes, %.3f seconds \n', elapsedTimeHours,elapsedTimeMinuts,elapsedTimeSeconds);
