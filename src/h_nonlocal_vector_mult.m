function  Hnlx = h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,X,S,kptvec)
% @brief   Calculates Hamiltonian vector product, where Vnl is the nonlocal
%          pseudopotential to be determined using the info stored in S.
%
% @authors  Abhiraj Sharma <asharma424@gatech.edu>
%           Qimen Xu <qimenxu@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% 
% @param DLii        Discrete laplacian component in 1D along ith direction
% @param DGi         Discrete gradient component in 1D along ith direction
% @param Veff        Effective potential N x 1 vector
% @param X           N x Ns matrix of eigenfunctions of Hamiltonian
% @param kptvec      k-point vector for the current Block diagonalized problem
% @param Hnlx        Hamiltonian times vector  
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%========================================================================================

% (-0.5*Lap + Veff) * X 
%Hnlx = -0.5*(lapVec(DL11,DL22,DL33,DG1,DG2,DG3,X,S)) + Veff * X;
Hnlx = -0.5*(lapVec(DL11,DL22,DL33,DG1,DG2,DG3,X,S)) + bsxfun(@times,Veff,X);

% Vnl * X
if (kptvec(1) == 0 && kptvec(2) == 0 && kptvec(3) == 0)
	fac = 1.0;
else
	fac = 1.0i;
end

for J = 1:S.n_atm
	Chi_X_mult = zeros(S.Atom(J).angnum,size(X,2));
	for img = 1:S.Atom(J).n_image_rc
		% Chi_X_mult = Chi_X_mult + (S.Atom(J).rcImage(img).Chi_mat .* S.W(S.Atom(J).rcImage(img).rc_pos))' * X(S.Atom(J).rcImage(img).rc_pos,:)* ...
		Chi_X_mult = Chi_X_mult + ( bsxfun(@times, S.Atom(J).rcImage(img).Chi_mat, S.W(S.Atom(J).rcImage(img).rc_pos)) )' * X(S.Atom(J).rcImage(img).rc_pos,:) * ...
					(exp(dot(kptvec,(S.Atoms(J,:)-S.Atom(J).rcImage(img).coordinates)*fac)));
	end
	% Chi_X_mult = Chi_X_mult .* S.Atom(J).gamma_Jl;
	Chi_X_mult = bsxfun(@times,Chi_X_mult, S.Atom(J).gamma_Jl);
	for img = 1:S.Atom(J).n_image_rc
		Chi = S.Atom(J).rcImage(img).Chi_mat * (exp(-dot(kptvec,(S.Atoms(J,:)-S.Atom(J).rcImage(img).coordinates)*fac)));
		Hnlx(S.Atom(J).rcImage(img).rc_pos,:) = Hnlx(S.Atom(J).rcImage(img).rc_pos,:) + Chi * Chi_X_mult;
	end
end
