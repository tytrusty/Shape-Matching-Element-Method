function test_normals
clf;
parts=nurbs_from_iges('rounded_cube.iges');
parts=nurbs_plot(parts);

N = nurbs_normals(parts{1}.srf.nurbs, parts{1}.UV, parts{1}.p);

% Create line segment along normal
p1 = parts{1}.x0;
p2 = parts{1}.x0 + 10*N';

% Plotting normals
plot3([p1(1,:); p2(1,:)],[p1(2,:); p2(2,:)],[p1(3,:); p2(3,:)], ...
    'LineWidth', 3);
end
