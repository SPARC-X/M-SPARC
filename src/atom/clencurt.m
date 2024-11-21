function [w] = clencurt(N)
% Clenshaw-Curtis Quadrature
% weights w for Clenshaw-Curtis quadrature for chebyshev grids in [-1,1]
% Adapted from "Spectral methods in MATLAB" by Llyod N Trefethen, SIAM (2000)
%===============================================================================
theta = pi*(0:N)'/N;
w = zeros(1,N+1);
ii = 2:N;
v = ones(N-1,1);
if mod(N,2)==0
    w(1) = 1/(N^2-1);
    w(N+1) = w(1);
    for k = 1:N/2-1
        v = v - 2*cos(2*k*theta(ii))/(4*k^2-1);
    end
    v = v - cos(N*theta(ii))/(N^2-1);
else
    w(1) = 1/N^2;
    w(N+1) = w(1);
    for k = 1:(N-1)/2
        v = v- 2*cos(2*k*theta(ii))/(4*k^2-1);
    end
end
w(ii) = 2*v/N;
end

