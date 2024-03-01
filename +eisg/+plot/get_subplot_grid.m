function grid = get_subplot_grid( n_plots )

factors = factor( n_plots );
factors = flip( factors );
n_factors = length( factors );
prod1_end_ind = floor( sqrt( n_factors ) );

n_cols = prod( factors(1:prod1_end_ind) ); % factors with odd index
n_rows = prod( factors(prod1_end_ind+1:end) ); % factors with even index
grid = [n_rows, n_cols];

end