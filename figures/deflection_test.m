function deflection_test

    
% Read in cube mesh 
files={'b1.igs', 'b2.igs', 'b3.igs', 'b4.igs', 'b5.igs', 'b6.igs','b7.igs','b8.igs','b9.igs','b10.igs'};%, ...
      % 'b7.igs', 'b8.igs'};%, 'b9.igs', 'b10.igs'};%, 'b20.igs'};
n = numel(files);
dofs=zeros(n,1);
nparts=zeros(n,1);
def=zeros(n,1);

YM = 1e6; %in Pascals
pr =  0.45;
[lambda, mu] = emu_to_lame(YM, pr);

options.dt = 0.1;
options.order = 2;
options.pin_function = @pin_function;
options.gravity = -10;
options.lambda = lambda;
options.mu = mu;
options.rho = 1e3;
options.distance_cutoff = 0.501;
% options.distance_cutoff = 1;
options.k_stability = 1e7;
options.enable_secondary_rays = 1;
options.fitting_mode = 'hierarchical';
options.save_output = 0;
options.save_obj = 0;
options.plot_com=1;
% options.plot_points = 1;
options.x_samples = 100;
options.y_samples = 3;
options.z_samples = 3;

iters=20;
for i=1:n
    parts = nurbs_from_iges(files{i});
    dofs(i) = numel(parts{1}.p)*numel(parts);
    [def(i),deflections] = deflection_test_step(parts,options);
    nparts(i) = numel(parts);
    figure(3);
    plot(1:iters, deflections);
    hold on;
end
figure(4);
plot(nparts, abs(def ./ -4.1579),'LineWidth',2);
xlabel('DOFs');
ylabel('m');
title('Max Deflection');
% yline(-4.1579,'Color','r');
end