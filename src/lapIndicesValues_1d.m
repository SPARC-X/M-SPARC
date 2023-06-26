function S = lapIndicesValues_1d(S)
% @ brief    Calculates laplacian and gradient (in 1D) indices' values without Bloch factor
% @ authors
%         Abhiraj Sharma <asharma424@gatech.edu>
%         Qimen Xu <qimenxu@gatech.edu>
%         Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%============================================================================
if S.cell_typ < 3
	Nx = S.Nx; Ny = S.Ny; Nz = S.Nz;
	n0 = S.FDn;
	w1 = S.w1;
	w2 = S.w2;
	dx = S.dx;
	dy = S.dy;
	dz = S.dz;

	% D_xx laplacian in 1D
	%-----------------------

	% Initial number of non-zeros: including ghost nodes
	nnzCount = (2 * n0 + 1) * Nx;

	% Row and column indices and the corresponding non-zero values
	% used to generate sparse matrix DL11 s.t. DL11(I(k),II(k)) = V(k)
	I = zeros(nnzCount,1);
	V = zeros(nnzCount,1);
	II = zeros(nnzCount,1);
	rowCount = 1;
	count = 1;
	coef_dxx = 1/dx^2;

	% Find non-zero entries that use forward difference
	for ii = 1:Nx
		% diagonal element
		I(count) = rowCount; II(count) = ii;
		V(count) = w2(1)*coef_dxx ;
		count = count + 1;
		% off-diagonal elements
		for q = 1:n0
			% ii + q
			I(count) = rowCount; II(count) = ii+q;
			V(count) = w2(q+1)*coef_dxx;
			count = count + 1;
			% ii - q
			I(count) = rowCount; II(count) = ii-q;
			V(count) = w2(q+1)*coef_dxx;
			count = count + 1;
			
		end
		rowCount = rowCount + 1;
	end

	if S.BCx == 1
		% Removing outside domain entries (for periodic code this is unnecessary)
		isIn = (II >= 1) & (II <= Nx);
		S.I_11 = I(isIn); S.II_11 = II(isIn); S.V_11 = V(isIn);
	elseif S.BCx == 0
		S.isOutl_11 = (II<1); S.isOutr_11 = (II>Nx); % Warning: Assumed influence of only neighboring cells
		S.I_11 = I; S.II_11 = mod(II+(Nx-1),Nx)+1; S.V_11 = V;
	end

	% D_yy laplacian in 1D
	%-----------------------

	% Initial number of non-zeros: including ghost nodes
	nnzCount = (2 * n0 + 1) * Ny;

	% Row and column indices and the corresponding non-zero values
	% used to generate sparse matrix DL22 s.t. DL22(I(k),II(k)) = V(k)
	I = zeros(nnzCount,1);
	V = zeros(nnzCount,1);
	II = zeros(nnzCount,1);
	rowCount = 1;
	count = 1;
	coef_dyy = 1/dy^2;

	% Find non-zero entries that use forward difference
	for ii = 1:Ny
		% diagonal element
		I(count) = rowCount; II(count) = ii;
		V(count) = w2(1)*coef_dyy;
		count = count + 1;
		% off-diagonal elements
		for q = 1:n0
			% ii + q
			I(count) = rowCount; II(count) = ii+q;
			V(count) = w2(q+1)*coef_dyy;
			count = count + 1;
			% ii - q
			I(count) = rowCount; II(count) = ii-q;
			V(count) = w2(q+1)*coef_dyy;
			count = count + 1;
			
		end
		rowCount = rowCount + 1;
	end

	if S.BCy == 1
		% Removing outside domain entries (for periodic code this is unnecessary)
		isIn = (II >= 1) & (II <= Ny);
		S.I_22 = I(isIn); S.II_22 = II(isIn); S.V_22 = V(isIn);
	elseif S.BCy == 0
		S.isOutl_22 = (II<1); S.isOutr_22 = (II>Ny); % Warning: Assumed influence of only neighboring cells
		S.I_22 = I;  S.II_22 = mod(II+(Ny-1),Ny)+1; S.V_22 = V;
	end

	% D_zz laplacian in 1D
	%-----------------------

	% Initial number of non-zeros: including ghost nodes
	nnzCount = (2 * n0 + 1) * Nz;

	% Row and column indices and the corresponding non-zero values
	% used to generate sparse matrix DL33 s.t. DL33(I(k),II(k)) = V(k)
	I = zeros(nnzCount,1);
	V = zeros(nnzCount,1);
	II = zeros(nnzCount,1);
	rowCount = 1;
	count = 1;
	coef_dzz = 1/dz^2;

	% Find non-zero entries that use forward difference
	for ii = 1:Nz
		% diagonal element
		I(count) = rowCount; II(count) = ii;
		V(count) = w2(1)*coef_dzz ;
		count = count + 1;
		% off-diagonal elements
		for q = 1:n0
			% ii + q
			I(count) = rowCount; II(count) = ii+q;
			V(count) = w2(q+1)*coef_dzz;
			count = count + 1;
			% ii - q
			I(count) = rowCount; II(count) = ii-q;
			V(count) = w2(q+1)*coef_dzz;
			count = count + 1;
			
		end
		rowCount = rowCount + 1;
	end

	if S.BCz == 1
		% Removing outside domain entries (for periodic code this is unnecessary)
		isIn = (II >= 1) & (II <= Nz);
		S.I_33 = I(isIn); S.II_33 = II(isIn); S.V_33 = V(isIn);
	elseif S.BCz == 0
		S.isOutl_33 = (II<1); S.isOutr_33 = (II>Nz); % Warning: Assumed influence of only neighboring cells
		S.I_33 = I; S.II_33 = mod(II+(Nz-1),Nz)+1; S.V_33 = V;
	end

	if S.cell_typ == 2
		% Create 1D gradient in all directions
		%---------------------------------------

		% x-direction
		%-------------
		nnz_x = 2*n0*Nx;
		G = zeros(nnz_x,1);
		R = zeros(nnz_x,1);
		A = zeros(nnz_x,1);
		rowCount = 1;
		count = 1;
		coef_dx = 1/dx;

		for ii = 1:Nx
			for q = 1:n0
				% ii + q
				G(count) = rowCount; R(count) = ii+q;
				A(count) = w1(q+1)*coef_dx;
				count = count + 1;
				% ii - q
				G(count) = rowCount; R(count) = ii-q;
				A(count) = -w1(q+1)*coef_dx;
				count = count + 1;
			end
			rowCount = rowCount + 1;
		end

		if S.BCx == 1
			% Removing outside domain entries (for periodic code this is unnecessary)
			isIn = (R >= 1) & (R <= Nx);
			S.I_1 = G(isIn); S.II_1 = R(isIn); S.V_1 = A(isIn);
		elseif S.BCx == 0
			S.isOutl_1 = (R<1); S.isOutr_1 = (R>Nx); % Warning: Assumed influence of only neighboring cells
			S.I_1 = G; S.II_1 = mod(R+(Nx-1),Nx)+1; S.V_1 = A;
		end

		% y-direction
		%-------------

		nnz_y = 2*n0*Ny;
		G = zeros(nnz_y,1);
		R = zeros(nnz_y,1);
		A = zeros(nnz_y,1);
		count =1;
		rowCount =1;
		coef_dy = 1/dy;

		for jj = 1:Ny
			for q = 1:n0
				% jj + q
				G(count) = rowCount; R(count) = jj+q;
				A(count) = w1(q+1)*coef_dy;
				count = count + 1;
				% jj - q
				G(count) = rowCount; R(count) = jj-q;
				A(count) = -w1(q+1)*coef_dy;
				count = count + 1;
			end
			rowCount = rowCount + 1;
		end

		if S.BCy == 1
			% Removing outside domain entries (for periodic code this is unnecessary)
			isIn = (R >= 1) & (R <= Ny);
			S.I_2 = G(isIn); S.II_2 = R(isIn); S.V_2 = A(isIn);
		elseif S.BCy == 0
			S.isOutl_2 = (R<1); S.isOutr_2 = (R>Ny);
			S.I_2 = G; S.II_2 = mod(R+(Ny-1),Ny)+1; S.V_2 = A;
		end


		% z-direction
		%-------------

		nnz_z = 2*n0*Nz;
		G = zeros(nnz_z,1);
		R = zeros(nnz_z,1);
		A = zeros(nnz_z,1);
		count =1;
		rowCount =1;
		coef_dz = 1/dz;

		for kk = 1:Nz
			for q = 1:n0
				% kk + q
				G(count) = rowCount; R(count) = kk+q;
				A(count) = w1(q+1)*coef_dz;
				count = count + 1;
				% kk - q
				G(count) = rowCount; R(count) = kk-q;
				A(count) = -w1(q+1)*coef_dz;
				count = count + 1;
			end
			rowCount = rowCount + 1;
		end

		if S.BCz == 1
			% Removing outside domain entries (for periodic code this is unnecessary)
			isIn = (R >= 1) & (R <= Nz);
			S.I_3 = G(isIn); S.II_3 = R(isIn); S.V_3 = A(isIn);
		elseif S.BCz == 0
			S.isOutl_3 = (R<1); S.isOutr_3 = (R>Nz);
			S.I_3 = G; S.II_3 = mod(R+(Nz-1),Nz)+1; S.V_3 = A;
		end
	end

elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
   S = lapIndicesValues_1d_cychel(S);

end