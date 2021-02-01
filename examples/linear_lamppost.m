function linear_lamppost

function pinned_ids = pin_function(x)
    verts_to_pin = 64; 
    [~,I] = mink(x(3,:),verts_to_pin);
    pinned_ids = I(1:3:verts_to_pin);
end

iges_file = 'lamppost.iges';

% To avoid singular nurbs jacobian with excessive pinning, I up the sample
% density to ensure we have enough unpinned samples.
sample_density=2;
part = nurbs_from_iges(iges_file, sample_density);

YM = 1e6; %in Pascals
pr = 0.4;
% YM = 2e11; %in Pascals
% pr = 0.32;
[lambda, mu] = emu_to_lame(YM, pr);

options.dt = 0.02;
options.order = 2;
options.gravity = 0;
options.rho = 3e3;
% options.rho = 1.27e3;
% options.rho = 8e3;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;
options.distance_cutoff = 1.0;
options.k_stability = 5e5;
options.plot_points = 1;
options.x_samples = 2;
options.y_samples = 5;
options.z_samples = 35;
options.save_iges=1;
cutoff_distance=cutoff_heuristic(part, 0.8);
options.distance_cutoff = 0.8;

options.f_external = [-1e1 0 0];
options.f_external_time = 2.0;

% vem_simulate_nurbs_newtons(part, options);
vem_simulate_nurbs(part, options);

end