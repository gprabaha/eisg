function labels = cluster2label(clusters,numElements)

numClusters = length(clusters);
labels = zeros(1,numElements);
for i = 1:numClusters
    labels(clusters{i}(:)) = i;
end