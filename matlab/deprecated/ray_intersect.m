function [dist, faces] = ray_intersect(xrange, coord, verts, face_bins, Octree)
    % Find all bins the x-axis aligned ray intersects with
    intersect_bins = all((coord > Octree.BinBoundaries(:,2:3) & ...
                          coord < Octree.BinBoundaries(:,5:6)),2);

    face_intersect = face_bins(intersect_bins);
    face_intersect = face_intersect(~cellfun(@isempty, face_intersect));
    
    % Intersect along all valid faces
    face_intersect = cell2mat(face_intersect);
    if isempty(face_intersect)
       dist = [];
       faces = [];
       return
    end
    P0=verts(face_intersect(:,1),:);
    P1=verts(face_intersect(:,2),:);
    P2=verts(face_intersect(:,3),:);
    or = [xrange(2) coord];
    
    D = [xrange(1) coord] - or;
    D = D ./ norm(D);
    [dist,fids] = ray_tri(P0, P1, P2, or, D);

    % Output faces if desired
    if nargout > 1
        faces = face_intersect(fids,:);
    end
    
    % Plotting intersection bins             
    %     boxH = OT.plot('Visible',0);
    %     for ii = 1:numel(valid_bins)
    %        i = valid_bins(ii);
    %        set(boxH(i),'Color','r','LineWidth', 1,'Visible',1)
    %        hold on;
    %     end
end