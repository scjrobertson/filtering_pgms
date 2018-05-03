function C = simplifiedMultiprod(A, B, d1, d2, d3, m)
% SIMPLIFIEDMULTIPROD - Calculates C = A*B for a 3D array.
%   C = simplifiedMultiprod(A, B, d1, d2, d3, m)
%
%   A simplified version of multiprod for 3D array equations of the form
%   C = A*B.
%
%   See also multiprod, rightMultiprod, leftMultiprod, quadraticMultiprod.
%
%   Inputs
%       A - (d1, d2, m) array. 
%       B - (d2, d3, m) array. 
%       d1 - scalar. The first dimension of the left matrix.
%       d2 - scalar. The second dimension of the left matrix.
%       d3 - scalar. The second dimension of the right matrix.
%       m - scalar. The number of matrices.
%
%   Outputs
%       C - (d1, d3, m) array. C = A*B.
D = reshape(A, [d1 d2 1 m]);
E = reshape(B, [1 d2 d3 m]);
F = sum(bsxfun(@times, D, E), 2);
C = reshape(F, [d1 d3 m]);
end