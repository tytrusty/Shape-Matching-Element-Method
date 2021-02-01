function  test_trimmed_nurbs

% Works ok, just need to make this more oop
% make all surfaces enclosed in a boundary surface...
% support lines in addition to curve boundaries

clf;
% data=iges2matlab('cylinder2.iges');
data=iges2matlab('wrench.iges');

subd=10;

function [points, normals] = process_boundaries(boundary_srf)
    % Mark surface as trimmed
    data{boundary_srf.sptr}.is_trimmed = 1;
    
    points = [];
    normals = [];
    
    for i=1:boundary_srf.n
        boundary = data{boundary_srf.bdpt(i)};
        for j=1:numel(boundary.pscpt)
            bound = data{boundary.pscpt{j}};
            if bound.type == 110 % Line
                subd = 5;
            elseif bound.type == 126 % NURBS CURVE
                subd = 50;
            end
            tt=linspace(bound.v(1), bound.v(2), subd);
            [p, dp] = nrbdeval(bound.nurbs, bound.dnurbs, tt); 
            plot(p(1,:), p(2,:),'-','LineWidth',2);
            hold on;


            unit_z = repmat([0 0 1]',1, subd);
            N = cross(dp,unit_z);
            N = N ./ vecnorm(N);
            points = [points p];
            normals = [normals N];
%             end
        end
    end
    
end

for ii=1:numel(data)  
   if data{ii}.type == 128
       data{ii}.is_trimmed = 0;
   end
end

for ii=1:numel(data)   
    if data{ii}.type == 143
        clf;
    	srf_boundary = data{ii};
        [points, N] = process_boundaries(srf_boundary);
        p2 = points + N*.1;
        plot([points(1,:); p2(1,:)], [points(2,:); p2(2,:)], 'Color', 'blue');
        
        % Sample a bunch of points on the surface
        srf_res = 21;
        srf = data{srf_boundary.sptr};
        u = linspace(srf.u(1), srf.u(2), srf_res);
        v = linspace(srf.v(1), srf.v(2), srf_res);
        [U,V] = meshgrid(u,v);
        UV = [U(:) V(:)];
        is_valid = ones(size(UV,1),1);
        T = delaunay(UV(:,1), UV(:,2));

        % Computing min distance from each point to the curve points
        diff_X = UV(:,1) - points(1,:);
        diff_Y = UV(:,2) - points(2,:);
        dist_sqr = diff_X.^2 + diff_Y.^2;
        
        % Computing if inside polygon via winding numbers.
        % This is pretty doo doo, I should just shoot a ray a check line
        % intersections along ray.
        winding = diff_X .* N(1,:) + diff_Y .* N(2,:); % numerator term
        winding = winding ./ (dist_sqr * 2 * pi);      % denominator
        winding = sum(winding, 2);
        e = 0.5;
        sgn = (winding+e)./e;

        is_valid = sgn < 0;
        
%         patch('Faces',T,'Vertices',UV);
        
        vids = find(is_valid);
        [toreplace, bywhat] = ismember(T, vids);
        toremove = ~all(toreplace, 2);
        T(toreplace) = bywhat(toreplace);
        T(toremove,:) = [];
        
        
        UV_valid = UV(is_valid,:);
        plot(UV(:,1), UV(:,2), '.', 'Color', 'r', 'MarkerSize', 10);
        hold on;
        plot(UV_valid(:,1), UV_valid(:,2), '.', 'Color', 'g', 'MarkerSize', 20);
        hold on
        patch('Faces',T,'Vertices',UV_valid);
        
%         p3 = UV_valid + min_point_dir(is_valid,:);
%         plot([UV_valid(:,1)'; p3(:,1)'], [UV_valid(:,2)'; p3(:,2)']);
        
        disp('Press a button to look at the next thing');
        waitforbuttonpress;
    end
end

end

