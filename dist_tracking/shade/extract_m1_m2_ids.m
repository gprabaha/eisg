function [m1, m2] = extract_m1_m2_ids(sessions)

tmp_f = fcat.from( sessions, {'session'} );
bfw.add_monk_labels( tmp_f );

m1 = cellstr( tmp_f, 'id_m1' );
m2 = cellstr( tmp_f, 'id_m2' );

end