function tf = above_sem(psth, n_devs)

sigma = plotlabeled.sem( psth' )';
mu = mean( psth, 2 );
tf = false( size(psth) );

for i = 1:size(psth, 1)
  tf(i, :) = psth(i, :) > mu(i, :) + sigma(i, :) .* n_devs;
end

end