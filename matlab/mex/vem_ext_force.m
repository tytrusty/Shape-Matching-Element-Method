% VEM_EXT_FORCE Compute external forces for a specificed constant force
% vector
%
% Inputs:
%   ext_f   3 by 1 force vector
%   mass    M by 1 list of per point masses
%   Y       M by 3kn gradient of deformed positions with respect to their
%                    polynomial coefficients.
%
% Outputs:
%    g  3kn by 1 force vector for the set of coefficients where k is the
%                size of the monomial basis and n is the number of parts.
%