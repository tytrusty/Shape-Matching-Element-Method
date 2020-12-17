function w = blending_weights(parts, X)

parts=nurbs_from_iges('rocket_4.iges',8,0);
parts=plot_nurbs(parts);

X = raycast_quadrature(parts, [8 8], 5);

% Triangulate nurbs patch
faces=[];
verts=[];
for ii=1:numel(parts)
    fvc = surf2patch(parts{ii}.plt,'triangles');
    fvc.faces = fvc.faces + size(verts,1);
    faces=[faces; fvc.faces];
    verts=[verts; fvc.vertices];   
end
clf;
FV.faces = faces;
FV.vertices = verts;
points = X;
[distances,surface_points] = point2trimesh(FV, 'QueryPoints', points, 'UseSubSurface', false);

diff = surface_points - points;
rays = - diff ./ vecnorm(diff, 2, 2);
X2 = X + rays * 10;

P0 = FV.vertices(FV.faces(:,1),:);
P1 = FV.vertices(FV.faces(:,2),:);
P2 = FV.vertices(FV.faces(:,3),:);

for i=1:size(X,1)
    if abs(distances(i)) < 1e-5
        %plot3(X(i,1),X(i,2),X(i,3),'.','Color','m','MarkerSize',20);
    else
        [dist,~] = ray_tri(P0,P1,P2, X(i,:), rays(i,:));
        dist(dist < 0) = [];
        dist=min(dist);
        X2 = X(i,:) + rays(i,:) * dist(1);
        if i==42
            plot3([X(i,1); X2(1)], [X(i,2); X2(2)], [X(i,3); X2(3)],'Color','r', 'LineWidth',3);
        end

    end
    hold on;

end


patch(FV,'FaceAlpha',.5); xlabel('x'); ylabel('y'); zlabel('z'); axis equal; hold on
plot3M = @(XYZ,varargin) plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),varargin{:});
plot3M(points(42,:),'*r')
plot3M(surface_points(42,:),'*k')
plot3M(reshape([shiftdim(points(42,:),-1);shiftdim(surface_points(42,:),-1);shiftdim(points(42,:),-1)*NaN],[],3),'g','LineWidth',3)


end