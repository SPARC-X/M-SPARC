function DG = blochGradient(S,kptvec,grad_dir)
% @ brief   Calculates gradient of a function along "grad_dir" direction
% @ authors
%         Abhiraj Sharma <asharma424@gatech.edu>
%         Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param kptvec          k-point vector for the current Block diagonalized problem
% @param grad_dir        direction along which gradient has to be calculated
% @param DG              Discrete gradient along grad_dir direction 
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%===================================================================================

if S.cell_typ < 3

	if (grad_dir == 1)
		I = S.G_I_1;
		J = S.G_J_1;
		V = S.G_V_1;
		BC = S.BCx;
		if BC == 0
			isOutl = S.G_JOutl_1;
			isOutr = S.G_JOutr_1;
			L = S.L1;
			kpt = kptvec(1);
		end
	elseif (grad_dir == 2)
		I = S.G_I_2;
		J = S.G_J_2;
		V = S.G_V_2;
		BC = S.BCy;
		if BC == 0
			isOutl = S.G_JOutl_2;
			isOutr = S.G_JOutr_2;
			L = S.L2;
			kpt = kptvec(2);
		end
	else
		I = S.G_I_3;
		J = S.G_J_3;
		V = S.G_V_3;
		BC = S.BCz;
		if BC == 0
			isOutl = S.G_JOutl_3;
			isOutr = S.G_JOutr_3;
			L = S.L3;
			kpt = kptvec(3);
		end
	end

	% Map the elements of discrete laplacian
	if BC == 0
		% Phase factor
		if (kpt == 0)
			phasefac_l = 1.0;
			phasefac_r = 1.0;
		else
			phasefac_l = exp(-1i*kpt*L);
			phasefac_r = exp(1i*kpt*L);
		end

		V(isOutl) =  V(isOutl) * phasefac_l;
		V(isOutr) =  V(isOutr) * phasefac_r;
	end

	DG = sparse(I,J,V,S.N,S.N);

else
	DG = blochGradient_cychel(S,kptvec,grad_dir);
end
