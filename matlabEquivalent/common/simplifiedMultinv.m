function B = simplifiedMultinv(A, d, m)
% SIMPLIFIEDMULTINV -- A simplified multinv for a very specific 3D array
% format.
%   B = simplifiedMultinv(A, d, m)
%
%   A simplified multidimensional array for a very specific 3D array
%   format.
%
%   See also multinv.
%
%   Inputs
%       A - (d, d, m) array. A 3D array of covariance matrices.
%       d - scalar. The dimension of the covariance matrices.
%       m - scalar. The number of covariance matrices.
%
%   Output
%       B - (d, d, m) array. A 3D array of precision matrices.
I = reshape(1:d*m, [d 1 m]);
I = repmat(I, [1 d 1]);
J = permute(I, [2 1 3]);
A = sparse(I(:), J(:), A(:));
RHS = repmat(eye(d), [m 1]);
B = reshape(A \ RHS, [d m d]);
B = permute(B, [1 3 2]);
end