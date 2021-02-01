function quadratic_rounded_cube

function pinned_ids = pin_function(x)
    % pinning corner points of the cube
    x_min = 0.4;
    y_min = 0.4;
    pinned_ids = find(x(1,:) > x_min & x(2,:) > y_min);
    
%         verts_to_pin = 9; 
%     [~,I] = mink(x(3,:),verts_to_pin);
%     pinned_ids = I(1:verts_to_pin);
end

iges_file = 'rounded_cube.iges';
part = nurbs_from_iges(iges_file);
YM = 1e5; %in Pascals
% YM = 1e6; %in Pascals
pr =  0.45;
[lambda, mu] = emu_to_lame(YM, pr);
options.order = 2;
options.rho = 1e2; % for compression example
options.rho = 1e1;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;

options.distance_cutoff=cutoff_heuristic(part, 0.9);

% uncomment to effectively perform single shape matching
% options.enable_secondary_rays = 0;
% options.distance_cutoff = 100;
% options.gravity = -400; % for compression

% options.save_obj=1;
vem_simulate_nurbs_cube_anim(part, options);

end