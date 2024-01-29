function [betas, ps] = extract_mdl_beta_coeffs(mdl)

%   There are two possible terms in the model, x1, x2
%   For each one, look for the coefficient by name in the RowNames
%   property of the coefficient table. If it isn't there, use NaN as the
%   value for the coefficient, otherwise extract it.

betas = nan( 1, 2 );
ps = nan( 1, 2 );

coeffs = mdl.Coefficients;
row_names = coeffs.Properties.RowNames;

%   @TODO: Handle the case that x1 is not in the table of coefficients, in
%   which case `x1_ind` will be empty ([])

%if strcmp(row_names, 'x1') = [];
 %  x1_ind = 'NaN';
x1_ind = find( strcmp(row_names, 'x1') );
if (~isempty(x1_ind)) 
    betas(1) = coeffs.Estimate(x1_ind);
    ps(1) = coeffs.pValue(x1_ind);
end 

x2_ind = find(strcmp(row_names,'x2') );
if (~isempty(x2_ind))
    betas(2) = coeffs.Estimate(x2_ind);
    ps(2) = coeffs.pValue(x2_ind);
end 

end