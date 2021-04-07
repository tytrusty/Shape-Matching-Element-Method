function scalability_test

fig=figure(2);
    
% Read in cube mesh 
iges_file = 'cube.iges';
parts = nurbs_from_iges(iges_file);


YM = 1e6; %in Pascals
pr =  0.45;
[lambda, mu] = emu_to_lame(YM, pr);


n=50;
nnz=zeros(n,1);
dofs=zeros(n,1);
time1=zeros(n,1);
time2=zeros(n,1);
substeps=1;
for i=1:n
    clf;
    parts=nurbs_plot(parts);
    options.x_samples = 5*i;
    options.y_samples = 3;
    options.z_samples = 3;
    options.gravity = -10;
    options.lambda = lambda;
    options.mu = mu;
    options.rho = 1e3;
    options.distance_cutoff = 0.501;

    options.order = 2; % shouldn't matter, right?
    options.distance_cutoff = 0.51;
    dofs(i) = numel(parts{1}.p)*numel(parts);
    for j=1:substeps
        [nnz(i), t1, t2] = scalability_test_step(parts,options);
        time1(i) = time1(i) + t1/substeps;
        time2(i) = time2(i) + t2/substeps;
    end
    drawnow

%     fn=sprintf('output/img/scalability_%03d.png',i);
%     saveas(fig,fn);
            
    % Push end farther, and duplicate 4 parts
    parts{1}.p(1:3:end) = parts{1}.p(1:3:end) + 1;
    parts{1}.x0 = squeeze(sum(parts{1}.J .* parts{1}.p,1));
    parts{1}.hires_x0 = squeeze(sum(parts{1}.hires_J .* parts{1}.p,1));

    for j=3:6 % the indices into the cube to 
        new_part=parts{j};
        new_part.p(1:3:end) = new_part.p(1:3:end) + i;
        new_part.x0 = squeeze(sum(new_part.J .* new_part.p,1));
        new_part.hires_x0 = squeeze(sum(new_part.hires_J .* new_part.p,1));
        parts{numel(parts)+1} = new_part;
    end
end
figure(3);
plot(dofs, nnz,'LineWidth',2);
xlabel('DOFs');
ylabel('%');
title('Percent nnz');
figure(4);
% plot(1:n, times,'LineWidth',2);
% area(dofs, [time1 time2],'LineWidth',2, 'EdgeAlpha',0.0);
plot(dofs, time1,'LineWidth',2);
xlabel('DOFs');
ylabel('time (s)');
title('Total time in seconds');
figure(5);
plot(dofs, time2,'LineWidth',2);
xlabel('DOFs');
ylabel('time (s)');
title('Total time in seconds');
% figure(6)
% area(dofs, [time1 time2],'LineWidth',2, 'EdgeAlpha',0.0);
% xlabel('DOFs');
% ylabel('time (s)');
% title('Total time in seconds');
end