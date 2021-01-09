clear;
cut_off = 100; % threshold
[d1, d2] = meshgrid(0:5:150, 0:5:150); 
dist_w = @(d1, d2) max(1.0 - d1 ./ min(d2, cut_off), 0.0);

surf(d1, d2, dist_w(d1, d2));
% colormap summer
alpha(0.5)
% shading interp
xlabel('Primary Ray Length')
ylabel('Secondary Ray Length')
zlabel('Distance Weight')