function surfaces=nurbs_sample(data, enable_trimming, sample_density)

nsurfaces = 0;

% Initially mark all surfaces as untrimmed
for ii=1:numel(data)  
	if data{ii}.type == 128
        data{ii}.is_trimmed = 0;
       
        s = sample_spline(data{ii}.s, data{ii}.m1, sample_density);
        t = sample_spline(data{ii}.t, data{ii}.m2, sample_density);
        [U,V] = meshgrid(s,t);
        UV = [U(:) V(:)]';
        data{ii}.UV = UV;
            
        nsurfaces = nsurfaces + 1;
	end
end

surfaces = cell(nsurfaces,1);
ii = 1;

if enable_trimming
    process_trimmed();
end
process_untrimmed();

function process_trimmed()
    % Samples points on all the trimmed surfaces
    for i=1:numel(data)   
        if data{i}.type == 143 && data{i}.boundarytype == 1
            [p1,p2,N] = boundary_lines(data{i});

            % Sample a bunch of points on the surface
            srf_ii = data{i}.sptr;

            % Mark surface as trimmed
            data{srf_ii}.is_trimmed = 1;
            UV = sample_rays(p1, p2, data{srf_ii});

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
            uv_struct.UV = data{i}.UV;
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
                nlines = 100;
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

function u_vals = sample_intersections(u_spans, nurbs, alpha)
    eps = 1e-12;
    degree = nurbs.m1;
    knots = unique(nurbs.s(degree+1:end-degree));
    interval = @(x) find(x-knots > -eps, 1, 'last');
    u_vals = [];

    for i = 1:size(u_spans,2)

        idx1 = interval(u_spans(1,i));
        idx2 = interval(u_spans(2,i));

        % Special case when both lie in the same span.
        if idx1==idx2
            samples = linspace(u_spans(1,i), u_spans(2,i), degree + alpha);
            u_vals = [u_vals samples];
        end

        for j = idx1:idx2-1
            % Sample along the knot span. Remove final samples if not at
            % end index so we don't have duplicates
            u1 = knots(j);
            u2 = knots(j+1);
            if j == idx1
                u1 = u_spans(1,i);
            elseif j == idx2-1
                u2 = u_spans(2,i);
            end
            samples = linspace(u1, u2, degree + alpha);
            if j < idx2 - 1
               samples = samples(1:end-1) ;
            end
            u_vals = [u_vals samples];
        end
    end
end

function UV = sample_rays(p1, p2, nurbs)
    UV = [];
    alpha = 3;

    origin_u = nurbs.u(1);
    v_range = [min([p1(2,:) p2(2,:)]) max([p1(2,:) p2(2,:)])];
    v_vals = sample_spline(nurbs.t, nurbs.m2, alpha);

    eps=1e-5;
    v_vals(1) = v_range(1)+eps;
    v_vals(end) = v_range(2)-eps;

%     clf;
%     plot([p1(1,:); p2(1,:)],[p1(2,:); p2(2,:)], 'LineWidth',5, 'Color',[0.5 0.5 0.5]);
%     hold on;

    for i = 1:numel(v_vals)
        % Find intersection points.
        isect = (p1(2,:) <= v_vals(i) & p2(2,:) >= v_vals(i)) ...
              | (p1(2,:) >= v_vals(i) & p2(2,:) <= v_vals(i));

        if nnz(isect) > 0
            origin = [origin_u v_vals(i)]';
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

            t_pnts = [];

            % If odd number of intersections, perform point sampling,
            % otherwise sample along the ray taking consideration of the
            % knots in the u direction.
            if rem(numel(t),2) == 1
                % Point sampling
                new_t_pnts = origin + t .* [1 0]';
                t_pnts = [t_pnts new_t_pnts];
            else
                u_spans = origin_u + t;
                u_spans = reshape(u_spans,2,[]);
                u_samples = sample_intersections(u_spans, nurbs, alpha);
                v_samples = repelem(v_vals(i), numel(u_samples));
                t_pnts = [t_pnts [u_samples; v_samples]];
            end
            UV = [UV t_pnts];

%             plot([nurbs.u(1) nurbs.u(2)], [v_vals(i) v_vals(i)],'-','Color','b','LineWidth',2);
%             hold on;
%             if numel(t_pnts) > 0
%                 plot(t_pnts(1,:), t_pnts(2,:), '.', 'Color', 'g', 'MarkerSize', 20);
%             end
        end
    end
end
end