function linear_lamppost

% % this will blow up somehow
% function pinned_ids = pin_function(x)
%     pinned_ids = find(x(3,:) < 0.1);
% end

function pinned_ids = pin_function(x)
    verts_to_pin = 20; 
    [~,I] = mink(x(3,:),verts_to_pin);
    pinned_ids = I(1:verts_to_pin);
end

iges_file = 'lamppost.iges';

part = nurbs_from_iges(iges_file);

YM = 1e5; %in Pascals
pr = 0.4;
% YM = 2e11; %in Pascals
% pr = 0.32;
[lambda, mu] = emu_to_lame(YM, pr);

% options.dt = 0.1;
options.order = 1;
options.gravity = 0;
options.rho = 1.27e3;
options.pin_function = @pin_function;
options.lambda = lambda;
options.mu = mu;
options.distance_cutoff = 1.0;
options.k_stability = YM * 1e2;
options.plot_points = 1;
options.x_samples = 10;
% options.y_samples = 2;
options.z_samples = 15;
options.f_external = [-1e1 0 0];

% vem_simulate_nurbs_newtons(part, options);
vem_simulate_nurbs(part, options);

end