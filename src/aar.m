function [x] = aar(A,b,x_guess,tol,max_iter,omega,beta,m,p,L,U)
% @brief	Preconditioned Alternating Anderson-Richardson (AAR) MATLAB 
%           code. Solves the system Ax = b. Last Modified: 26 March 2018 
%
% @param A       : square matrix (size N x N),
% @param b       : right hand side column vector (size N x 1),
% @param x_guess : initial guess (of size b),
% @param tol     : convergence tolerance
% @param max_iter: maximum number of iterations
% @param omega   : relaxation parameter (for Richardson update)
% @param beta    : extrapolation parameter (for Anderson update)
% @param m       : Anderson history, no. of previous iterations to be considered in extrapolation
% @param p       : Perform Anderson extrapolation at every p th iteration
% @param L       : Lower part of LU Preconditioner
% @param L       : Upper part of LU Preconditioner
% @output x      : solution vector
%
% NOTE  : A and b are to be provided. Other input parameters can be passed as "[]", to use defaults.
% e.g.  : Run code as, tic;x=aar(A,b,[],[],[],[],[],[],[],[]);toc;
%
% @ref  : Suryanarayana, P., Pratapa, P. P., & Pask, J. E. (2019). 
%         Alternating Anderson–Richardson method: An efficient alternative 
%         to preconditioned Krylov methods for large, sparse linear systems. 
%         Computer Physics Communications, 234, 278-285.
%
% @authors	Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%			Phanisri Pradeep Pratapa 
%           (In collaboration with John E. Pask.)
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech

 
N=length(b);
DX=zeros(N,m);
DF=zeros(N,m);
nb = 1/norm(b) ;
if isempty(x_guess)
	x_guess = ones(N,1);  % initial guess vector
end
x_prev = x_guess;
relres = 1+tol;
count = 1;
while count<=max_iter && relres > tol
	res1 = (b-A*x_prev) ;
	relres = norm(res1)*nb;
	res = U\(L\(res1)); 
	
	% STORE HISTORY
	if count>1
		DX(:,mod(count-2,m)+1) = x_prev-Xold;
		DF(:,mod(count-2,m)+1) = res-Fold;
	end
	Xold = x_prev;
	Fold = res;
	
	% UPDATE ITERATE
	if rem(count,p)~=0   % RICHARDSON UPDATE  
		x_new = x_prev + omega*res; 
	else % ANDERSON UPDATE, apply every p iters
		x_new = x_prev + beta*res - (DX + beta*DF)*(pinv(DF'*DF)*(DF'*res));
	end
	
	x = x_prev;
	x_prev = x_new;
	count = count + 1;
end

if count-1 == max_iter
	fprintf(' AAR exceeded maximum iterations and converged to a relative residual of %g. \n',relres);
else
	fprintf(' AAR converged to a relative residual of %g in %d iterations.\n',relres,count-1);
end

end
