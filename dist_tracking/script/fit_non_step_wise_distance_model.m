function mdl = fit_non_step_wise_distance_model(X, y, offset)

mdl = fitglm(...
    X, y, 'linear', 'Link','log','Distribution' ,'poisson','Offset',log(offset) );

end