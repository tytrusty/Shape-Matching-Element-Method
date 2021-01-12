clear;
cut_off = 50; % threshold
[d1, total] = meshgrid(0:2:150, 0:2:150); 
dist_w_func = @(d1, d2) (max(1.0 - d1 ./ min(d2, cut_off), 0.0)) .* (d1 <= d2);
dist_w = dist_w_func(d1, total-d1);
% dist_w(d1 > total - d1) = NaN;
% dist_w(total - d1 < 0) = NaN;

t = surf(d1, total, dist_w);
colormap summer
alpha(0.5)
set(t, 'edgealpha', 0.6)
% shading interp
xlabel('Primary Ray Length')
ylabel('Total Ray Length')
zlabel('Distance Weight')