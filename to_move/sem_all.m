function s = sem_all(v)

s = std( v(:) ) / sqrt( numel(v) );

end