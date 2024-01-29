X = rand(10, 2);
y = randi(20, size(X, 1), 1);
offset = rand( size(X, 1), 1 ) * 0.2;

mdl1 = stepwiseglm(...
    X, y, 'constant','Lower','linear','Upper','linear' ...
    , 'Link','log','Distribution' ,'poisson','Offset',log(offset), 'Verbose', 0);

mdl2 = fitglm(...
    X, y, 'linear', 'Link','log','Distribution' ,'poisson','Offset',log(offset) );