function Y = chebyshev_filter(X,m,a,b,a0,DL11,DL22,DL33,DG1,DG2,DG3,Veff,S,kpt_vec)
% @brief    CHEBYSHEV_FILTER performs chebyshev fitering on the given
%           states.
%
% @param X  Current states to be filtered.
% @param m  Chebyshev polynomial degree.
% @param a  Lower bound of the filter.
% @param b  Upper bound of the filter.
% @param a0 Cutoff of the filter.
%
% @ref:  
% Y. Zhou, et al, 2016. Parallel Self-Consistent-Field Calculations via
% Chebyshev-Filtered Subspace Acceleration.
% https://www-users.cs.umn.edu/~saad/PDF/umsi-2006-101.pdf
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech


% Note: the nonlocal pseudopotential is taken into account
e = (b-a)/2;
c = (b+a)/2;
sigma = e/(a0 - c); sigma1 = sigma;
gamma = 2/sigma1;
% HX = H * X;
HX = h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,X,S,kpt_vec);
Y = (sigma1/e)*(HX-c*X);

ii = 2;
while(ii <= m)
	sigma2 = 1/(gamma - sigma);
	% HX = H * Y;
	HX = h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,Y,S,kpt_vec);
	Ynew = (2*sigma2/e)*(HX - c*Y) - ((sigma*sigma2)*X);
	X = Y;
	Y = Ynew;
	sigma = sigma2;
	ii = ii + 1;
end

end
