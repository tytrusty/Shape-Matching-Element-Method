function t_new = sample_spline(u, degree, alpha)
    if nargin < 3
        alpha = 2; % 1 sample in the range (midpoint)
    end
    u = u(degree+1:end-degree);
    t_new = [];
    for i=1:numel(u)-1
        s_b = u(i);
        s_e = u(i+1);
        samples = linspace(s_b, s_e, degree + alpha);
        t_new = [t_new samples(1:end-1)];
    end
    t_new = [t_new samples(end)];
end