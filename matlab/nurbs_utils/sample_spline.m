% function [t,w] = sample_spline(u, degree, alpha)
%     u = u(degree+1:end-degree);
%     t = [];
%     w = [];
%     for i=1:numel(u)-1
%         s_b = u(i);
%         s_e = u(i+1);
%         N = degree + alpha;
%         ti = linspace(s_b, s_e, N);
%         t = [t ti(1:end-1)];
%         wi = (s_e-s_b)/N;
%         w = [w repelem(wi, N-1)];
%     end
%     t = [t ti(end)];
%     w = [w wi];
%     fprintf('numel %d ... sum %f\n', numel(w), sum(w));
% end
function [t,w] = sample_spline(u, degree, alpha)
    u = u(degree+1:end-degree);
    t = [];
    w = [];
    for i=1:numel(u)-1
        s_b = u(i);
        s_e = u(i+1);
        N = degree + alpha;
        [samples,w_new]=lglnodes(N, s_e, s_b);
        t = [t; samples];
        w = [w; -w_new];
    end
    t = t';
    w = w';
    %fprintf('numel %d ... sum %f\n', numel(w), sum(w));
end