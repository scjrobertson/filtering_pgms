function C = leftMultiprod(A, B, d1, d2, d3, m)
% LEFTMULTIPROD - Calculates C = A*B for a 3D array.
%   C = leftMultiprod(A, B, d1, d2, d3, m)
%
%   A simplified version of multiprod for 3D array equations of the form
%   C = A*B.
%
%   See also multiprod, rightMultiprod, quadraticMultiprod.
%
%   Inputs
%       A - (d1, d2) array. 
%       B - (d2, d3, m) array. A 3D array of the targets' covariance matrices. 
%       d1 - scalar. The first dimension of the left matrix.
%       d2 - scalar. The second dimension of the left matrix.
%       d3 - scalar. The second dimension of the right matrix.
%       m - scalar. The number of matrices.
%
%   Outputs
%       C - (d1, d3, m) array. C = A*B.
F = reshape(B, [d2 d3*m]);
C = reshape(A*F, [d1 d3 m]);
end