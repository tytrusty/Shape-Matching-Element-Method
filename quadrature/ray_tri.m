function [dist_out,IDs_out] = ray_tri(P0, P1, P2, or, D)
% Code source: https://www.mathworks.com/matlabcentral/fileexchange/49160-fast-mesh-mesh-intersection-using-ray-tri-intersection-with-octree-spatial-partitioning
% Vectorized ray-triangle intersection algorithm of Muller and Trumbore 
% (1997)

% Author: Thomas Seers, the Univerity of Manchester, UK.

epsilon = 0.00001;
flagV = ones(size(P0,1),1);
E1 = P1-P0;
E2 = P2-P0;
% D=repmat(D,size(E2,1),1);
% Q = cross(D,E2,2);
Q = [(D(:,2).*E2(:,3)-D(:,3).*E2(:,2)) (D(:,3).*E2(:,1)-D(:,1).*E2(:,3)) (D(:,1).*E2(:,2)-D(:,2).*E2(:,1))];
A = sum(E1.*Q,2);

flagV(A > -epsilon & A < epsilon) = 0;

F = ones(size(P0,1),1) ./ A;
S = or-P0;
U = F.*sum(S.*Q,2);

flagV(U<0.0) = 0;

% R = cross(S,E1,2);
R = [(S(:,2).*E1(:,3)-S(:,3).*E1(:,2)) (S(:,3).*E1(:,1)-S(:,1).*E1(:,3)) (S(:,1).*E1(:,2)-S(:,2).*E1(:,1))];
V = F.*sum(D.*R,2);

%flagV(V<0.0 | U+V>1.0) = 0;
flagV(V<epsilon | U+V>(1.0-epsilon)) = 0;

T = F.*sum(E2.*R,2); 

% extract outputs and exclude the cloned results from triangles found in
% multiple bins;
IDs_out = find(flagV);
dist_out = T(IDs_out,:);

clear E1 E2 Q A F S U V T
end

