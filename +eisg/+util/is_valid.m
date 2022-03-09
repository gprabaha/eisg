function unit_validity = is_valid(sorted_neural_data, file_ind, unit_ind, params)

unit_validity = false;
if params.include_unsure_units
    validity_cond = sorted_neural_data(file_ind).validity(unit_ind)==1 | isnan(sorted_neural_data(file_ind).validity(unit_ind));
else
    validity_cond = sorted_neural_data(file_ind).validity(unit_ind)==1;
end
if validity_cond
    unit_validity = true;
end

end