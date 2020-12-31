function surfaces=nurbs_sample(data, enable_trimming)

nsurfaces = 0;

% Initially mark all surfaces as untrimmed
for ii=1:numel(data)  
	if data{ii}.type == 128
        data{ii}.is_trimmed = 0;
       
        s = sample_spline(data{ii}.s, data{ii}.m1);
        t = sample_spline(data{ii}.t, data{ii}.m2);
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

function t_new = sample_spline(u, degree)
    u = u(degree+1:end-degree);
    %u = unique(u);
    t_new = [];
    for i=1:numel(u)-1
        s_b = u(i);
        s_e = u(i+1);
        samples = linspace(s_b, s_e, degree + 2);
        t_new = [t_new samples(1:end-1)];
    end
    t_new = [t_new samples(end)];
end

function process_trimmed()
    % Samples points on all the trimmed surfaces
    for i=1:numel(data)   
        if data{i}.type == 143 && data{i}.boundarytype == 1
            [p1,p2,N] = boundary_lines(data{i});

            % Sample a bunch of points on the surface
            srf_ii = data{i}.sptr;

            % Mark surface as trimmed
            data{srf_ii}.is_trimmed = 1;

            dist = point_line_min_distances(data{srf_ii}.UV, p1, p2);
            [~,I] = min(dist, [], 2);

            % Check angles between nearest point ray and normal.
            diff_p = data{srf_ii}.UV - p1(:, I);
            angles = dot(diff_p, N(:,I));
            ToKeep = angles < 0;
%             clf;
%             plot(data{srf_ii}.UV(1,:), data{srf_ii}.UV(2,:),'.','Color','g');
%             hold on;
%             plot(data{srf_ii}.UV(1,ToRemove), data{srf_ii}.UV(2,ToRemove),'.','Color','r');
            UV = data{srf_ii}.UV(:,ToKeep);

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

function dist = point_line_min_distances(origins, p1, p2)
    % ref: http://paulbourke.net/geometry/pointlineplane/
    diff_X = origins(1,:)' - p1(1,:);
    diff_Y = origins(2,:)' - p1(2,:);
    dp = p2 - p1;
    norm_sqr = dot(dp,dp,1);

    t = (diff_X .* dp(1,:) + diff_Y .* dp(2,:)) ./ norm_sqr;
    t(:) = max(min(t(:), 1), 0);

    min_X = p1(1,:) + t .* dp(1,:);
    min_Y = p1(2,:) + t .* dp(2,:);

    diff_X = origins(1,:)' - min_X;
    diff_Y = origins(2,:)' - min_Y;
    dist = diff_X.^2 + diff_Y.^2;
end

end