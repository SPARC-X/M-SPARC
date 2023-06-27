function S = gradIndicesValues(S)
if S.cell_typ < 3
	Nx = S.Nx; Ny = S.Ny; Nz = S.Nz;
	N = S.N;
	n0 = S.FDn;
	w1 = S.w1;
	dx = S.dx;
	dy = S.dy;
	dz = S.dz;

	% Initial number of non-zeros: including ghost nodes
	nnz_count = 2*n0*N ;

	% Row numbers and non-zero values
	I = zeros(nnz_count,1) ;
	V = zeros(nnz_count,1) ;

	% Indices of the columns
	II = zeros(nnz_count,1);
	JJ = zeros(nnz_count,1) ;
	KK = zeros(nnz_count,1) ;

	% Gradient along x_direction
	row_count = 1;
	count = 1 ;
	for kk=1:Nz
		for jj=1:Ny
			for ii=1:Nx
				% off-diagonal elements
				for p=1:n0
					% ii+p
					I(count) = row_count; II(count) = ii+p; JJ(count) = jj; KK(count) = kk;
					V(count) = w1(p+1)/dx;
					count = count + 1;
					% ii-p
					I(count) = row_count; II(count) = ii-p ; JJ(count) = jj; KK(count) = kk;
					V(count) = -w1(p+1)/dx;
					count = count + 1;
				end
				row_count = row_count+1;
			end
		end
	end

	if S.BCx == 1
		% Removing outside domain entries (for periodic code this is unnecessary)
		isIn = (II >= 1) & (II <= Nx);
		I = I(isIn); II = II(isIn); JJ = JJ(isIn); KK = KK(isIn); V = V(isIn);
	elseif S.BCx == 0
		S.G_JOutl_1 = (II<1); S.G_JOutr_1 = (II>Nx); % Warning: Assumed influence of only neighboring cells
		II = mod(II+(Nx-1),Nx)+1;
	end

	% Getting linear indices of the columns
	S.G_J_1 = (KK-1)*Nx*Ny + (JJ-1)*Nx + II;
	S.G_I_1 = I;
	S.G_V_1 = V;

	% Gradient along y_direction
	row_count = 1;
	count = 1 ;
	for kk=1:Nz
		for jj=1:Ny
			for ii=1:Nx
				% off-diagonal elements
				for p=1:n0
					% ii+p
					I(count) = row_count; II(count) = ii; JJ(count) = jj+p; KK(count) = kk;
					V(count) = w1(p+1)/dy;
					count = count + 1;
					% ii-p
					I(count) = row_count; II(count) = ii ; JJ(count) = jj-p; KK(count) = kk;
					V(count) = -w1(p+1)/dy;
					count = count + 1;
				end
				row_count = row_count+1;
			end
		end
	end

	if S.BCy == 1
		% Removing outside domain entries (for periodic code this is unnecessary)
		isIn = (JJ >= 1) & (JJ <= Ny);
		I = I(isIn); II = II(isIn); JJ = JJ(isIn); KK = KK(isIn); V = V(isIn);
	elseif S.BCy == 0
		S.G_JOutl_2 = (JJ<1); S.G_JOutr_2 = (JJ>Ny); % Warning: Assumed influence of only neighboring cells
		JJ = mod(JJ+(Ny-1),Ny)+1;
	end

	% Getting linear indices of the columns
	S.G_J_2 = (KK-1)*Nx*Ny + (JJ-1)*Nx + II;
	S.G_I_2 = I;
	S.G_V_2 = V;

	% Gradient along z_direction
	row_count = 1;
	count = 1 ;
	for kk=1:Nz
		for jj=1:Ny
			for ii=1:Nx
				% off-diagonal elements
				for p=1:n0
					% ii+p
					I(count) = row_count; II(count) = ii; JJ(count) = jj; KK(count) = kk+p;
					V(count) = w1(p+1)/dz;
					count = count + 1;
					% ii-p
					I(count) = row_count; II(count) = ii ; JJ(count) = jj; KK(count) = kk-p;
					V(count) = -w1(p+1)/dz;
					count = count + 1;
				end
				row_count = row_count+1;
			end
		end
	end

	if S.BCz == 1
		% Removing outside domain entries (for periodic code this is unnecessary)
		isIn = (KK >= 1) & (KK <= Nz);
		I = I(isIn); II = II(isIn); JJ = JJ(isIn); KK = KK(isIn); V = V(isIn);
	elseif S.BCz == 0
		S.G_JOutl_3 = (KK<1); S.G_JOutr_3 = (KK>Nz); % Warning: Assumed influence of only neighboring cells
		KK = mod(KK+(Nz-1),Nz)+1;
	end

	% Getting linear indices of the columns
	S.G_J_3 = (KK-1)*Nx*Ny + (JJ-1)*Nx + II;
	S.G_I_3 = I;
	S.G_V_3 = V;

else
	S = gradIndicesValues_cychel(S);
end


