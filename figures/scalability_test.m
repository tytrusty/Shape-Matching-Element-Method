function scalability_test

fig=figure(2);
    
% Read in cube mesh 
iges_file = 'cube.iges';
parts = nurbs_from_iges(iges_file);

n=50;
nnz=zeros(n,1);
times=zeros(n,1);
for i=1:n
    clf;
    parts=nurbs_plot(parts);
    
    options.x_samples = 5*i;
    options.y_samples = 3;
    options.z_samples = 3;
    options.distance_cutoff = 0.51;

    [nnz(i), times(i)] =scalability_test_step(parts,options)
    drawnow
    fn=sprintf('output/img/scalability_%03d.png',i);
    saveas(fig,fn);
            
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
plot(1:n, nnz,'LineWidth',2);
xlabel('# cubes');
ylabel('%');
title('Percent nnz');
figure(4);
plot(1:n, times,'LineWidth',2);
xlabel('# cubes');
ylabel('time (s)');
title('Total time in seconds');
end