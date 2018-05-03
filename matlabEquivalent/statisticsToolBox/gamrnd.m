function sample = gamrnd(a,b,nrow,ncol)

% Purpose: 
% Generate random numbers from gamma distribution
% -----------------------------------
% Density:
% f(x) = 1/(G(a) * b^a) * x^(a-1) * exp(-x/b)
% E(X) = a*b, Var(X) = a*b^2
% -----------------------------------
% Algorithm: 
% Squeeze rejection sampling proposed by
% Marsaglia, G., and W. W. Tsang. 
% "A Simple Method for Generating Gamma Variables."
% -----------------------------------
% Usage:
% a = degree of freedom parameter (shape parameter) 
% b = scale parameter
% nrow = number of rows (if return a vector or matrix of random numbers)
% ncol = number of columns (if return a vector or matrix of random numbers)
% -----------------------------------
% Returns:
% sample = random numbers from Gamma(a,b)
% -----------------------------------
% Notes:
% Support vector or matrix input of parameter a,b
% In that case, no need to specify nrow and ncol.
% It will return a vector or matrix of random number with conformable size.
%
% Written by Hang Qian, Iowa State University
% Contact me:  matlabist@gmail.com

if nargin == 4
    a = a .* ones(nrow,ncol);
    b = b .* ones(nrow,ncol);
end

if nargin == 3
    a = a .* ones(nrow,1);
    b = b .* ones(nrow,1);
end

numel_a = numel(a);
numel_b = numel(b);
if numel_a >= numel_b    
    ndraws = numel_a;
    d_vec = a - 1/3;
    c_vec = 1./sqrt(9*d_vec);
    sample = a;
else
    ndraws = numel_b;
    d_vec = (a - 1/3) * ones(ndraws,1);
    c_vec = 1./sqrt(9*d_vec);
    sample = b;
end


for r = 1:ndraws    
    d = d_vec(r);
    c = c_vec(r);
    while 1
        x = randn;
        v = 1 + c*x;
        if v < 0
            continue
        end
        v = v^3;
        U = rand;
        if U < 1 - 0.0331 * x^4 || log(U) < 0.5*x^2 + d*(1 - v + log(v) )
            sample(r) = d * v;
            break
        end
    end
end

sample = sample .* b;
