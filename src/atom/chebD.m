function [D,r] = chebD(N,R)
% Chebyshev Differentiation
% Compute D = differentiation matrix, r = grids
% Modified from "Spectral methods in MATLAB" by Llyod N Trefethen, SIAM (2000)
% R is used to scale [-1,1] to [0,R]
%===============================================================================
Rmax = R;

% Chebyshev Grid
r = cos(pi*(0:N)/(N))';
r = Rmax*r + Rmax;

c = [2;ones(N-1,1);2].*(-1).^(0:N)';
X = repmat(r,1,N+1);
dX = X - X';
D = (c*(1./c)')./(dX + (eye(N+1)));
D = D - diag(sum(D'));
end

