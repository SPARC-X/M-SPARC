function  Hnlx = h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,X,S,kptvec,spin)
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

if S.usefock > 1
    Vexx = evaluateExactExchangePotential(S,X,kptvec,spin);
    Hnlx = Hnlx + S.hyb_mixing*Vexx;
end

if (S.xc == 4) && (S.countPotential > 0) % metaGGA, set a flag to seperate it from the 1st PBE SCF computation
    if S.nspin == 1
        VxcScan3 = S.VxcScan3;
    else % spin polarization mGSGA
        VxcScan3 = S.VxcScan3(:, spin); 
    end
    nCol = size(X, 2);
    if S.cell_typ == 2 % unorthogonal cell
        lapc_T = [S.lapc_T(1,1), S.lapc_T(2,1), S.lapc_T(3,1);
            S.lapc_T(2,1), S.lapc_T(2,2), S.lapc_T(3,2);
            S.lapc_T(3,1), S.lapc_T(3,2), S.lapc_T(3,3)];
        v3grad1 = blochGradient(S,kptvec,1) *X; 
        v3grad2 = blochGradient(S,kptvec,2) *X; 
        v3grad3 = blochGradient(S,kptvec,3) *X; 

        v3gradpsiTheKpt = [v3grad1(:), v3grad2(:), v3grad3(:)];
        v3gradpsiTheKptMLapT = v3gradpsiTheKpt*lapc_T;
        v3gradpsiTheKptMLapT_1 = VxcScan3 .* reshape(v3gradpsiTheKptMLapT(:, 1), S.N, nCol);
        v3gradpsiTheKptMLapT_2 = VxcScan3 .* reshape(v3gradpsiTheKptMLapT(:, 2), S.N, nCol);
        v3gradpsiTheKptMLapT_3 = VxcScan3 .* reshape(v3gradpsiTheKptMLapT(:, 3), S.N, nCol);
        Hnlx = Hnlx - 0.5*(blochGradient(S,kptvec,1)*v3gradpsiTheKptMLapT_1 + blochGradient(S,kptvec,2)*v3gradpsiTheKptMLapT_2 + blochGradient(S,kptvec,3)*v3gradpsiTheKptMLapT_3);
    else % orthogonal cell
        v3grad1 = VxcScan3 .* (blochGradient(S,kptvec,1) *X);
        v3grad2 = VxcScan3 .* (blochGradient(S,kptvec,2) *X);
        v3grad3 = VxcScan3 .* (blochGradient(S,kptvec,3) *X);
        Hnlx = Hnlx - 0.5*(blochGradient(S,kptvec,1)*v3grad1 + blochGradient(S,kptvec,2)*v3grad2 + blochGradient(S,kptvec,3)*v3grad3);
    end
end

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
