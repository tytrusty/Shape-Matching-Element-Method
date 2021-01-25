function quadratic_castle_cut_with_collision

function pinned_ids = pin_function(x)
    verts_to_pin = 3; 
    [~,I] = mink(x(3,:),verts_to_pin);
    pinned_ids = I(1:verts_to_pin);
%     pinned_ids = [];
end

iges_file = 'castle_simple_smaller_flip_cut.iges';

% To avoid singular nurbs jacobian with excessive pinning, I up the sample
% density to ensure we have enough unpinned samples.
sample_density=3;

part = nurbs_from_iges(iges_file, sample_density);

%material properties
YM = 1e4; %in Pascals
pr = 0.45;
[lambda, mu] = emu_to_lame(YM, pr);

options.dt = 0.005;
options.order = 2;
options.pin_function = @pin_function;
options.gravity = -10;
options.lambda = lambda;
options.mu = mu;
options.rho = 1.5e3;
options.distance_cutoff = 0.1;
options.k_stability = 1e5;
options.enable_secondary_rays = 1;
options.fitting_mode = 'hierarchical';
options.save_output = 0;
options.plot_points = 1;
options.x_samples = 12;
options.y_samples = 12;
options.z_samples = 15;

options.save_obj = true;
options.save_obj_path = 'output/obj_castle_cut_1e4_2balls/';

options.k_stability = YM*1e5;

% options.initial_velocity = [0 0 -10];
options.collision_with_other = false;
options.collision_other_position = [0 0 -10];
options.collision_with_other_sim = false;
options.self_collision = false;
options.collision_with_plane = true;
options.collision_plane_z = 0.684;
options.collision_with_sphere = true;
options.collision_ratio = 2.5e5;

options.collision_sphere_r = 0.0175;
options.collision_sphere_rho = 2e3;
options.collision_sphere_initial_veloity = [0 0 -10];

 options.collision_sphere_c = [0.3 0.5 0.8;
%                               0.3 0.4 0.8;
%                               0.4 0.5 0.8;
%                               0.4 0.3 0.8;
%                               0.4 0.4 0.8;
%                               0.4 0.5 0.8;
%                               0.5 0.3 0.8;
%                               0.5 0.4 0.8;
                              0.5 0.5 0.8];                           
options.collision_sphere_c = options.collision_sphere_c + 0.025 * rand(size(options.collision_sphere_c, 1), 3) ...
                              - 0.0125 - repmat([0 -0.05 0], size(options.collision_sphere_c, 1), 1);

                              
vem_simulate_nurbs_with_collision(part, options);

end