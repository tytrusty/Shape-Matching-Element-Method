function quadratic_castle_with_collision

function pinned_ids = pin_function(x)
    verts_to_pin = 10; 
    [~,I] = mink(x(3,:),verts_to_pin);
    pinned_ids = I(1:verts_to_pin);
%     pinned_ids = [];
end

% iges_file = 'lamppost.iges';
iges_file = 'castle_simple.iges';

% To avoid singular nurbs jacobian with excessive pinning, I up the sample
% density to ensure we have enough unpinned samples.
sample_density=1;

part = nurbs_from_iges(iges_file, sample_density);

%material properties
% YM = 5e3; %in Pascals
YM = 5e4; %in Pascals
pr = 0.45;
[lambda, mu] = emu_to_lame(YM, pr);

options.order = 2;
options.pin_function = @pin_function;
options.gravity = -10;
options.lambda = lambda;
options.mu = mu;
options.rho = 1.27e3;
% options.distance_cutoff = 1.5;
% options.distance_cutoff = 0.9;
options.distance_cutoff = 1;
options.k_stability = 1e5;
options.enable_secondary_rays = 1;
options.fitting_mode = 'hierarchical';
options.save_output = 0;
% options.plot_points = 1;
% options.x_samples = 10;
% options.y_samples = 10;
% options.z_samples = 10;
% vem_simulate_nurbs_newtons(part, options);

% options.initial_velocity = [0 0 -10];
options.collision_with_other = false;
options.collision_other_position = [0 0 -10];
options.collision_with_other_sim = false;
options.self_collision = false;
options.collision_with_plane = false;
options.collision_with_sphere = true;
options.collision_ratio = 5e4;

options.collision_sphere_r = 2.0;
options.collision_sphere_c = [-10 10 85;
                              -10 30 85;
                              -10 50 85;
                              -30 10 85;
                              -30 30 85;
                              -30 50 85;
                              -50 10 85;
                              -50 30 85;
                              -50 50 85];
                              

% options.collision_sphere_c = [0.2 0.2 1.5;
%                               0.2 -0.2 1.5;
%                               -0.2 0.2 1.5;
%                               -0.2 -0.2 1.5];

vem_simulate_nurbs_with_collision(part, options);

end