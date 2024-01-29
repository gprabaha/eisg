function center_roi_m2_eyes = get_hemifield_zero_point(offsets, roi_m2_eyes)

% Offsets_m1=offsets.m1;
Offsets_m2=offsets.m2;
Offsets_m2_m1=offsets.m2_to_m1;

if ( offsets.m1(1) == 1024) 
    addtl_offset = 0;
else
    addtl_offset = 1024;
end

% ROI info for eyes (m1 roi is m2's eyes from m1's perspective)
% roi_m1_eyes([1,3])=roi_m1_eyes([1,3])+Offsets_m1(1);
% roi_m1_eyes([2,4])=roi_m1_eyes([2,4])+Offsets_m1(2);
roi_m2_eyes([1,3])=Offsets_m2_m1(1)-(roi_m2_eyes([1,3])+Offsets_m2(1));
roi_m2_eyes([2,4])=roi_m2_eyes([2,4])+Offsets_m2(2);

% Get the center of ROI for eyes
% center_roi_m1_eyes=[(roi_m1_eyes(1)+roi_m1_eyes(3))/2 (roi_m1_eyes(2)+roi_m1_eyes(4))/2];
center_roi_m2_eyes=[(roi_m2_eyes(1)+roi_m2_eyes(3))/2 (roi_m2_eyes(2)+roi_m2_eyes(4))/2];

center_roi_m2_eyes(1) = center_roi_m2_eyes(1) - 1024 + addtl_offset;

end