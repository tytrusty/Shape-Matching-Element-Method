function surfaces=nurbs_sample_old(data, nrays, samples_per_ray, untrimmed_res)

nsurfaces = 0;

% Initially mark all surfaces as untrimmed
for ii=1:numel(data)  
   if data{ii}.type == 128
       data{ii}.is_trimmed = 0;
       nsurfaces = nsurfaces + 1;
   end
end

surfaces = cell(nsurfaces,1);
ii = 1;
process_trimmed();
process_untrimmed();
disp('done');

function process_trimmed()
    % Samples points on all the trimmed surfaces
    for i=1:numel(data)   
        if data{i}.type == 143 && data{i}.boundarytype == 1
            [p1,p2,N] = boundary_lines(data{i});

            % Sample a bunch of points on the surface
            srf_ii = data{i}.sptr;

            % Mark surface as trimmed
            data{srf_ii}.is_trimmed = 1;

            v_range = [min([p1(2,:) p2(2,:)]) max([p1(2,:) p2(2,:)])];
            v_vals = linspace(v_range(1)+1e-4, v_range(2)-1e-4, nrays);
            u_range = [data{srf_ii}.u(1) data{srf_ii}.u(2)];
            UV = sample_ray(p1, p2, v_vals, u_range, nrays, samples_per_ray);                      
            uv_struct.UV = UV;
            uv_struct.surf_ptr = srf_ii;
            uv_struct.line_1 = p1;
            uv_struct.line_2 = p2;
            uv_struct.line_N = N;
            uv_struct.is_trimmed = 1;
            surfaces{ii} = uv_struct;
            ii = ii + 1;
        end
    end
end

% Samples points on all the untrimmed surfaces
function process_untrimmed()
    for i=1:numel(data)   
        if data{i}.type == 128 && ~data{i}.is_trimmed
            srf = data{i};
            
            ratio = srf.ratio;
            res = max(round(untrimmed_res .*ratio ./ max(ratio)), 2);
            
            [s, ~, s_ic] = unique(data{i}.s);
            s_bins = histcounts(s_ic);
            s_new = [];
            for j=1:numel(s_bins)-1
                s_b = s(j);
                s_e = s(j+1);
                s_n = max(s_bins(j), s_bins(j+1));
                samples = linspace(s_b, s_e, s_n + 2);
                s_new = [s_new samples(1:end-1)];
            end
            s_new = [s_new samples(end)];
            s=s_new;
            
            
            [t, ~, t_ic] = unique(data{i}.t);
            t_bins = histcounts(t_ic);
            t_new = [];
            for j=1:numel(t_bins)-1
                t_b = t(j);
                t_e = t(j+1);
                t_n = max(t_bins(j), t_bins(j+1));
                samples = linspace(t_b, t_e, t_n + 2);
                t_new = [t_new samples(1:end-1)];
            end
            t_new = [t_new samples(end)];
            t=t_new;
            u = linspace(srf.u(1), srf.u(2), res(1));
            v = linspace(srf.v(1), srf.v(2), res(2));
%             [U,V] = meshgrid(u,v);
            [U,V] = meshgrid(s,t);
            UV = [U(:) V(:)]';
            uv_struct.UV = UV;
            uv_struct.surf_ptr = i;
            uv_struct.is_trimmed = 0;
            surfaces{ii} = uv_struct;
            ii = ii + 1;
        end
    end
end

function [p1,p2,N] = boundary_lines(boundary_srf)
    p1 = []; p2 = [];
    
    for i=1:boundary_srf.n
        boundary = data{boundary_srf.bdpt(i)};
        for j=1:numel(boundary.pscpt)
            curve = data{boundary.pscpt{j}};
            if curve.type == 110 % Line
                nlines = 1;
            elseif curve.type == 126 % NURBS CURVE
                nlines = 50;
            end
            tt = linspace(curve.v(1), curve.v(2), nlines+1);
            p = nrbeval(curve.nurbs, tt); 
            p1 = [p1 p(1:2,1:end-1)];
            p2 = [p2 p(1:2,2:end)];
        end
    end
    % Compute normals for each line segment sampled on boundary
    line = p2-p1;
    zero_z = repelem(0, size(p1,2));
    unit_z = repmat([0 0 1]',1, size(p1,2));
    N = cross([line; zero_z],unit_z);
    N = N ./ vecnorm(N);
    N=N(1:2,:);
end

function UV = sample_ray(p1, p2, v_vals, u_range, nrays, samples_per_ray)
    UV = [];
    for jj = 1:nrays
        % Find intersection points.
        isect = (p1(2,:) <= v_vals(jj) & p2(2,:) >= v_vals(jj)) ...
              | (p1(2,:) >= v_vals(jj) & p2(2,:) <= v_vals(jj));

        if nnz(isect) > 0
            origin = [u_range(1) v_vals(jj)]';
            p1_isect = p1(:,isect);
            p2_isect = p2(:,isect);

            % reference for ray-line intersection:
            % https://rootllama.wordpress.com/2014/06/20/ray-line-segment-intersection-test-in-2d/
            v1 = origin - p1_isect;
            v2 = p2_isect - p1_isect;
            v3 = repmat([0 1]', 1, size(v1,2));

            t = (v2(1,:).*v1(2,:) - v2(2,:).*v1(1,:)) ./ dot(v2,v3,1);
            t = sort(t);
            t = uniquetol(t);
            assert(rem(numel(t),2) == 0, ...
                'Number of intersection is not even! Trimming failed.')

            t = reshape(t,2,[]);

            t_pnts = [];
            for kk = 1:size(t,2)
                t_samples = linspace(t(1,kk), t(2,kk), samples_per_ray);
                new_t_pnts = origin + t_samples .* [1 0]';
                t_pnts = [t_pnts new_t_pnts];
            end
            UV = [UV t_pnts];
        end
    end
end

end