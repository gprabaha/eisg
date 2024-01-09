function data = standardization(data)
% data -- observations x dismensions
% standardization: every column is standardized with [0,1]

    nrow = size(data,1);
    colmin = min(data);
    colmax = max(data);
    dmax = colmax - colmin;
    data = data - repmat(colmin,nrow,1);
    data = data./repmat(dmax,nrow,1);