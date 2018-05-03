function C = rightMultiprod(B, A, d1, d2, d3, m)
% RIGHTMULTIPROD - Calculates C = B*A for a 3D array.
%   C = rightMultiprod(B, A, d1, d2, d3, m)
%
%   A simplified version of multiprod for 3D array equations of the form
%   C = B*A.
%
%   See also multiprod, leftMultiprod, quadraticMultiprod.
%
%   Inputs
%       B - (d1, d2, m) array. A 3D array of the targets' covariance matrices. 
%       A - (d2, d3) array. 
%       d1 - scalar. The first dimension of left matrices.
%       d2 - scalar. The second dimension of left matrics
%       d3 - scalar. The dimension of right matrices.
%       m - scalar. The number of matrices.
%
%   Outputs
%       C - (d1, d3, m) array. C = B*A.
D = reshape(permute(B, [1 3 2]), [d1*m d2]);
E = reshape(D*A, [d1 m d3]);
C = permute(E, [1 3 2]);
end