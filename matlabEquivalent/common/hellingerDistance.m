function distance = hellingerDistance(muX, Sxx,  muY, Syy)
% HELLINGERDISTANCE - The Hellinger distance between two Gaussians.
%   distance = hellingerDistance(muX, muY, SXX, SYY)
%
%   Calculate the Hellinger distance between to Gaussian distributions
%       
%   Input
%       mux - (n, 1) array. The 1st distrubtion's mean.
%       Sxx - (n, n) array. The 1st distrubtion's covariance.
%       muy - (n, 1) array. The 2nd distrubtion's mean.
%       Sxx - (n, n) array. The 2nd distrubtion's covariance.
%
%   Output
%       distance - scalar. The Hellinger distrubtion between two Gaussian
%           distrubtion.
muZ = muX - muY;
Szz = Sxx + Syy;

delta = sqrt(sqrt(det(Sxx)*det(Syy))/det(0.5*Szz));
epsilon = exp(-0.25*(muZ'/Szz)*muZ);

distance = sqrt(1 - delta*epsilon);
end