function [bins, face_bins, OT] = octree_mesh(faces, verts, bin_capacity)
    
    % Create octree structure
    if nargin < 3
        OT = OcTree(verts);
    else
        OT = OcTree(verts,'binCapacity',bin_capacity);
    end
    
    % Column vector of vertex indices 1:N
    ids = (1:size(verts,1))';

    % assign points to bins then collapse
    bins = cell(OT.BinCount,1);
    bin_ind = cell(OT.BinCount,1);
    for i = 1:numel(bins)
        bins{i} = verts(OT.PointBins==i,:);
        bin_ind{i} = ids(OT.PointBins==i,:);
    end
    
    % Place faces into octree bins
    face_bins = cell(size(bins,1),1);
    newIDBins = cell(size(bins,1),1);
    nfaces = size(faces,1);
    all_faces = vertcat(faces(:,1), faces(:,2), faces(:,3));
    for i=1:numel(bins)
        logi = ismember(all_faces, bin_ind{i});
        logicut = logi(1:nfaces) + logi(nfaces+1:nfaces*2) + logi((nfaces*2)+1:end);
        logicut(logicut > 0)=1;

        % Adds the vertex set of each face that intersects with this bin.
        face_bins{i} = faces(logical(logicut),:); 

        % For all faces that interesct with this bin, add all of their
        % associated vertices + vertex ids to the vertex bin and vertex id bin
        inpoints = unique(face_bins{i});
        inpoints(ismember(inpoints, bin_ind{i})) = [];     % must not be in set
        bins{i} = vertcat(bins{i}, verts(inpoints,:));     % add new vertices
        bin_ind{i} = vertcat(bin_ind{i}, ids(inpoints,:)); % add new vert ids

%         newIDBins{i} = (1:size(bin_ind{i},1))';
%         bin_size = size(face_bins{i},1);
%         face_stack = vertcat(face_bins{i}(:,1), face_bins{i}(:,2), face_bins{i}(:,3));
% 
%         % Remap global vertex ids to local ids within the single bin.
%         % These IDS in newIDBins serve as indices into bin_ind, which tell
%         % us the global vertex ids.
%         tempInds = arrayfun(@(x) find(bin_ind{i} == x,1,'first'), face_stack);
%         face_stack = newIDBins{i}(tempInds);
% 
%         % face_bins now contains local indices.
%         face_bins{i} = [face_stack(1:bin_size) ...
%                         face_stack(bin_size+1:bin_size*2)  ...
%                         face_stack((bin_size*2)+1:end)];
    end
end