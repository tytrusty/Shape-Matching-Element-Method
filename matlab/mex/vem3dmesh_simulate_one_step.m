% VEM2DMESH_NEOHOOKEAN_Q Compute the energy of Neohookean potential
% energy for a virtual element mesh
%
%
% Inputs:
%   A       2 by k shape matching matrix where k is the size of the
%                  monomial vectors
%   dF_dq   M by (2*2 * 2*N) gradient of deformation energy (F) with
%                            respect to the elements' configuration
%                            vectors (N is the # of points for this
%                            element)
%   dM_dX   M by (k*2) energy of the monomial basis vector with respect
%                      to the undeformed positions.
%   w       M by 1 list of weights
%   volume  M by 1 list of per point volumes
%   [C,D]   M by 2 material parameters matrix. 
%   k       scalar size of the monomial vectors
%   N       scalar the # of points for this element
%
% Outputs:
%    g  2*N by 1 energy vector
%
% Example:
%   % mesh (V,T)
%   vol = volume(V,T); %using gptoolbox volume command
%   dphidX = linear_tetmesh_dphi_dphidX(V,T);
%   q = reshape(V', 3*size(V,1), 1) %initial mesh position
%   YM = 2e6; %in Pascals
%   pr =  0.45
%   [lambda, mu] = emu_to_lame(YM*ones(size(T,1),1), pr*ones(size(T,1),1));
%   g = vem2dmesh_neohookean_dq(V,T,q,dphidX,vol,[0.5*mu,0.5*lambda]);