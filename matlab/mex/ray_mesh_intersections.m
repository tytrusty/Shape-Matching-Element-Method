% RAY_MESH_INTERSECTIONS  Find all hits (if it exists) for each ray.
%
% [flag, t] = ray_mesh_intersections(src, dir, V, F, k);
%
% Input:
%    src #rays by 3 list of 3D vector ray origins
%    dir #rays by 3 list of 3D vector ray directions
%    V  #V by 3 list of vertex positions
%    F  #F by 3 list of triangle indices
% Output:
%    id  #rays by k list of indices into F (0 if no hit)
%    t  #rays by k list of distances from the ray origin (inf if no hit)
%
