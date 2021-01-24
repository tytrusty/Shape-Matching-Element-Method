function quadratic_castle_with_collision

function pinned_ids = pin_function(x)
    verts_to_pin = 10; 
    [~,I] = mink(x(3,:),verts_to_pin);
    pinned_ids = I(1:verts_to_pin);
%     pinned_ids = [];
end

% iges_file = 'lamppost.iges';
iges_file = 'castle_simple_smaller_flip.iges';

% To avoid singular nurbs jacobian with excessive pinning, I up the sample
% density to ensure we have enough unpinned samples.
sample_density=1;

part = nurbs_from_iges(iges_file, sample_density);

%material properties
% YM = 5e3; %in Pascals
YM = 1e5; %in Pascals
pr = 0.45;
[lambda, mu] = emu_to_lame(YM, pr);

options.dt = 0.01;
options.order = 2;
options.pin_function = @pin_function;
options.gravity = -10;
options.lambda = lambda;
options.mu = mu;
options.rho = 2e3;
% options.distance_cutoff = 1.5;
% options.distance_cutoff = 0.9;
options.distance_cutoff = 0.1;
options.k_stability = 1e5;
options.enable_secondary_rays = 1;
options.fitting_mode = 'hierarchical';
options.save_output = 0;
% options.plot_points = 1;
% options.x_samples = 10;
% options.y_samples = 10;
% options.z_samples = 10;
% vem_simulate_nurbs_newtons(part, options);

options.save_obj = true;
options.save_obj_path = 'output/obj_castle/';

options.k_stability = YM*1e5;

% options.initial_velocity = [0 0 -10];
options.collision_with_other = false;
options.collision_other_position = [0 0 -10];
options.collision_with_other_sim = false;
options.self_collision = false;
options.collision_with_plane = true;
options.collision_plane_z = 0.684;
options.collision_with_sphere = true;
options.collision_ratio = 1.2e5;

options.collision_sphere_r = 0.0175;
options.collision_sphere_rho = 2e3;
options.collision_sphere_initial_veloity = [0 0 -10];
% options.collision_sphere_initial_veloity = [-15 0 -10];
% options.collision_sphere_c = [-1 -1 8.5;
%                               -1  0 8.5;
%                               -1  1 8.5;
%                                0 -1 8.5;
%                                0  0 8.5;
%                                0  1 8.5;
%                                1 -1 8.5;
%                                1  0 8.5;
%                                1  1 8.5];
% options.collision_sphere_c = [-2 -2 2;
%                               -2  0 2;
%                               -2  2 2;
%                                0 -2 2;
%                                0  0 2;
%                                0  2 2;
%                                2 -2 2;
%                                2  0 2;
%                                2  2 2];
options.collision_sphere_c = [0.3 0.3 0.8;
                              0.3 0.4 0.8;
                              0.3 0.5 0.8;
                              0.4 0.3 0.8;
                              0.4 0.4 0.8;
                              0.4 0.5 0.8;
                              0.5 0.3 0.8;
                              0.5 0.4 0.8;
                              0.5 0.5 0.8];
options.collision_sphere_c = options.collision_sphere_c + 0.025 * rand(9, 3) ...
                              - 0.0125 - repmat([0 -0.05 0], 9, 1);
% options.collision_sphere_c = options.collision_sphere_c + repmat([0.5 0 0], 9, 1);
                              

% options.collision_sphere_c = [0.0.5 0.0.5 0.5;
%                               0.0.5 -0.0.5 0.5;
%                               -0.0.5 0.0.5 0.5;
%                               -0.0.5 -0.0.5 0.5];

vem_simulate_nurbs_with_collision(part, options);

end