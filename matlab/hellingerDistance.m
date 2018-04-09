function distance = hellingerDistance(muX, muY, SXX, SYY)
    % hellinger_distance - The Hellinger distance between two Gaussians.
    muZ = muX - muY; 
    SZZ = SXX + SYY;
    
    delta = sqrt(sqrt(det(SXX)*det(SYY))/det(0.5*SZZ));
    epsilon = exp(-0.25*(muZ'/SZZ)*muZ);
    
    distance = sqrt(1 - delta*epsilon);
end