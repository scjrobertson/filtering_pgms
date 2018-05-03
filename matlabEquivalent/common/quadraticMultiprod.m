function C = quadraticMultiprod(B, A, Atranspose, d, m)
% QUADRATICMULTIPROD - Calculates C = A*B*A' for a 3D array.
%   C = quadraticMultiprod(B, A, Atranspose, d, m)
%
%   A simplified version of multiprod for 3D array equations of the form
%   C = A*B*A'. A very specific function for the prediction of multiple
%   targets covariance matrices.
%
%   Inputs
%       B - (d, d, m) array. A 3D array of the targets' covariance
%           matrices. d is the dimension of the state space.
%       A - (d, d) array. The state transition matrix.
%       Atranspose - (d, d) array. The transpose of the state transition
%           matrix.
%       d - scalar. The state space dimension.
%       m - scalar. The current number of targets.
%
%   Outputs
%       C - (d, d, m) array. C = A*B*A'.
D = reshape(permute(B, [1 3 2]), [d*m d]);
E = reshape(D*Atranspose, [d m d]);
F = reshape(permute(E, [1 3 2]), [d d*m]);
C = reshape(A*F, [d d m]);
end