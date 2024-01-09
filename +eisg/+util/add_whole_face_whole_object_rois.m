function events = add_whole_face_whole_object_rois(events)

[events.labels, ind] = bfw.make_whole_face_roi( events.labels );
events.events = events.events(ind, :);
[events.labels, ind] = bfw.make_whole_right_nonsocial_object_roi( events.labels );
events.events = events.events(ind, :);

end