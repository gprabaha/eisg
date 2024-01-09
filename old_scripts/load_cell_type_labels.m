function l = load_cell_type_labels(data_p)

error('!! Old function, check new implementation !!');

if ( nargin < 1 )
  data_p = '/Users/prabaha/repositories/eisg/processed_data';
end

l = load( fullfile(data_p, 'celltype_labels_p2v_combined.mat') );
l = l.(char(fieldnames(l)));

end