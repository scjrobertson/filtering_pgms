addpath common;

numberOfSimulations = 10000;
elapsedTime = zeros(3, numberOfSimulations);

d = 4;

model = generateModel(20, 0.95);
A = model.A;
Atranspose = model.Atranspose;
B = Atranspose*A;
R = model.R;

for i = 1:numberOfSimulations
    clear S;
    m = randi([1 1000]);
    S = rand(d, d, m);
    K = zeros(d, d, m);
    
    loopBeginTime = tic;
    for j = 1:m
        K(:, :, j) = A*S(:, :, j)*Atranspose + R;
    end
    elapsedTime(1, i) = toc(loopBeginTime);
    
    multiprodBeginTime = tic;
    foo = multiprod(S, Atranspose);
    Z = multiprod(A, foo) + R;
    elapsedTime(2, i) = toc(multiprodBeginTime);
    
    simplifiedBeginTime = tic;
    C = reshape(permute(S, [1 3 2]), [d*m d]);
    M = reshape(C*Atranspose, [d m d]);
    F = reshape(permute(M, [1 3 2]), [d d*m]);
    N = reshape(A*F, [d d m]) + R; 
    elapsedTime(3, i) = toc(simplifiedBeginTime);
end


meanTime = mean(elapsedTime, 2);

