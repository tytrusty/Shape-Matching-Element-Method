function linear_rounded_cube

function pinned_ids = pin_function(x)
    verts_to_pin = 4; 
    [~,I] = mink(x(3,:),verts_to_pin);
    pinned_ids = I(1:verts_to_pin);
end

iges_file = 'rounded_cube.iges';

% Resolution indicates how many point samples we will take on each
% e.g. 6 means we have 6 samples in both the U & V coordinates, so
%      a total of 36 samples across the NURBs patch.
resolution = 6;
part = nurbs_from_iges(iges_file, resolution, 0);

options.order = 1;
options.pin_function = @pin_function;

vem_simulate_nurbs(part, options);

end