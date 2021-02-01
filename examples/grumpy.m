function not_gumby

function pinned_ids = pin_function(x)
    %gumby_2
%     pinned_ids = find(x(1,:) >  -0.05 & x(1,:) < 0.4 & x(3,:) > 1.8 & x(3,:) < 3 & x(2,:) > 0);
    % gumby 4
%     pinned_ids = find(x(1,:) >  -.3 & x(1,:) < 0.4 & x(3,:) > 1.5 & x(3,:) < 3.3 & x(2,:) > 0);
%     excl = find(x(1,:) < 0 & x(3,:) < 1.6);

    % For the arm
    pinned_ids = find(x(1,:) >  -.3 & x(1,:) < 0.4 & x(3,:) > 1.5 & x(3,:) < 3.3);
    excl = find(x(3,:) < 2.35 & x(3,:) > 2.1 & x(1,:) > 0.1223);
    
%     pinned_ids(pinned_ids & excl) = [];
    pinned_ids = setdiff(pinned_ids, excl);
    pinned_ids = pinned_ids(1:3:end);
% pinned_ids=[];
end

iges_file = 'grumpy.iges';

% iges_file = 'rounded_cube.iges'; 

part = nurbs_from_iges(iges_file);
% YM = 1e4; %in Pascals
YM =  1e4; % leg 1e3; %in Pascals
pr =  0.4;

cutoff_distance=cutoff_heuristic(part, 0.75);

[lambda, mu] = emu_to_lame(YM, pr);
options.order = 2;
options.rho = 1e3;
options.gravity=0;
options.distance_cutoff = 1;%cutoff_distance;
options.distance_cutoff = 2;%cutoff_distance;
options.distance_cutoff = 0.6;%cutoff_distance;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;
options.plot_points = 1;
options.plot_com = 1;
options.k_stability=YM;
options.k_stability=1e9;
options.k_stability= 1e2; % leg 1e2;
options.enable_secondary_rays=1;
options.collision_ratio = 1e3;
% options.y_samples = 1;
% options.x_samples = 2;
% options.z_samples = 25;
options.y_samples = 1;
options.x_samples = 2;
options.z_samples = 41;
options.save_obj=1;
vem_sim_gumby(part, options);

end