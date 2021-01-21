function linear_trimmed_mug

function pinned_ids = pin_function(x)
    verts_to_pin = 4; 
    [~,I] = maxk(x(1,:),verts_to_pin);
    pinned_ids = I(1:verts_to_pin);
end

iges_file = 'mug_complex.iges';

part = nurbs_from_iges(iges_file);

YM = 2e3; %in Pascals
pr =  0.45;
[lambda, mu] = emu_to_lame(YM, pr);

options.rho = 1e-5;
options.order = 1;
options.pin_function = @pin_function;
options.gravity = -10;
options.lambda = lambda;
options.mu = mu;
% options.plot_points = 1;
options.fitting_mode = 'hierarchical';

options.x_samples = 3;
options.y_samples = 5;
options.z_samples = 5;

options.distance_cutoff = 100;
options.dt = 0.1;

vem_simulate_nurbs_newtons(part, options);
% vem_simulate_nurbs(part, options);
  
end