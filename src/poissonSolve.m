function [S] = poissonSolve(S, poisson_tol, Isguess)
% @brief    POISSONSOLVE solves the poisson equation for the 
%           electrostatic potential.
%
% @param poisson_tol    Tolerance for solving the poisson equation
%                       using iterative method AAR. 
% @param Isguess        1: guess vector provided,
%                       0: no guess vector available.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%

t1 = tic;

if S.cell_typ < 3
	f = poisson_RHS(S);
else 
	f = -4*pi*(S.rho(:,1)+S.b);
end

%fprintf(' Time taken for poisson_RHS: %f s\n',toc(t1));

% solve the poisson equation using a linear solver (AAR or gmres)
% AAR reference: see AAR.m 
if Isguess
	phi_guess = S.phi;
else
	phi_guess = [];
end

%[S.phi,conv_flag, relres, iter] = gmres(S.Lap_std,f,50,poisson_tol,50,S.LapPreconL,S.LapPreconU,phi_guess);
S.phi = aar(S.Lap_std,f,phi_guess,poisson_tol,S.MAXIT_POISSON,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU);

if(S.BC == 2)
	% To make sure integral phi is 0 (this removes the arbitariness of the
	% constant in phi calculation)
	S.phi = S.phi - dot(S.W,S.phi)/sum(S.W);
end
 
% [phi,conv_flag] = pcg(S.Lap,f,S.poisson_tol,50,S.LapPreconL,S.LapPreconU);
% fprintf('\n Poisson problem took %fs and %d iterations\n',toc(t1),iter(1)*iter(2));
% assert(conv_flag==0,'Electrostatic Poisson solve did not converge')
% S_Debug.relax(S.Relax_iter).poisson_flag_init = conv_flag;

fprintf(' Poisson problem took %fs\n',toc(t1));

end


function f = poisson_RHS(S)
% @brief	Poisson_RHS evaluates the right hand side of the poisson equation, 
%           including the boundary condtions, i.e. f = -4 * pi * ( rho + b - d) 
%           for cluster system, while for periodic system it's just 
%           f = -4 * pi * (rho + b).
rho = S.rho(:,1);
f = -4 * pi * (rho + S.b);

% for charged systems, add a uniform background charge so that total charge
% is 0
if S.NetCharge ~= 0
	unif_bkgd_chrg = S.NetCharge / (S.L1*S.L2*S.L3);
	f = f + (-4 * pi) * unif_bkgd_chrg;
	% fprintf(2,'total charge + background charge = %d\n',dot(S.W,rho+S.b+unif_bkgd_chrg));
end

if(S.BC == 1)
	% For cluster systems, we need to include boundary conditions d
	% RR = S.RR_AUG(S.isIn);
	for l = 0:S.l_cut
		multipole_moment(l+1).Qlm = sum(repmat(S.RR.^l .* (rho + S.b) .* S.W,1,2*l+1).* S.SH(l+1).Ylm )';
	end

	% Calculate phi using multipole expansion
	phi = zeros(size(S.RR_AUG_3D));
	for l = 0 : S.l_cut
		denom = (2*l+1)*S.RR_AUG_3D.^(l+1);
		for m = -l : l
			Ylm_AUG_3D = reshape(S.SH(l+1).Ylm_AUG(:,m+l+1),size(phi));
			phi = phi + Ylm_AUG_3D .* multipole_moment(l+1).Qlm(m+l+1) ./ denom;
		end
	end
	phi = 4 * pi * phi;
	phi(S.isIn) = 0;
elseif S.BC == 2
	return;
elseif S.BC == 3
	% find boundary conditions for periodic in 2D, dirichlet in 1D
	% Calculate phi
    phi = zeros(S.Nx+2*S.FDn, S.Ny+2*S.FDn, S.Nz+2*S.FDn);

	cellsize  = [S.L1, S.L2, S.L3];
	gridsizes = [S.Nx, S.Ny, S.Nz];
	meshsizes = [S.dx, S.dy, S.dz];
	bcs       = [S.BCx,S.BCy,S.BCz];

	% find which direction has Dirichlet BC
	dir_Z = find(bcs == 1); 
	dir_X = mod(dir_Z, 3) + 1;
	dir_Y = mod(dir_X, 3) + 1;

	% once we find the direction, we assume that direction is the Z
	% direction, the other two directions are then called X, Y
	NX = gridsizes(dir_X); NY = gridsizes(dir_Y); NZ = gridsizes(dir_Z);
	LX = cellsize(dir_X);  LY = cellsize(dir_Y);  LZ = cellsize(dir_Z);
	dX = meshsizes(dir_X); dY = meshsizes(dir_Y); dZ = meshsizes(dir_Z);
	A_XY = LX * LY; % area of (x',y') surface

	% reshape rho to 3D 
	rho = reshape(rho+S.b, S.Nx, S.Ny, S.Nz); % note here after rho = rho + b

	% permute rho so that the new z' direction has Dirichlet BC
	new_order = [dir_X, dir_Y, dir_Z]; % a permutation of [1,2,3]
	[~, reverse_order] = sort(new_order); % reverse order to get back
	rho = permute(rho,new_order);
	phi = permute(phi,new_order);
	% rho = permute(rho,reverse_order); % this will recover original

	% find fft of rho in the new X, Y directions
	% rho is a 3D matrix, this applies 2D fft to each dim higher than 2
	rho_hat = fft2(rho); 
	rho_hat(1,1,:) = 0; % remove the 0 frequency term (constant term in real-space)

	% rho_zp_bar = zeros(Nzp,1);
	% sum over x' and y' directions, \int (\rho) dxdy
	%rho_zp_bar = sum(sum(rho, dir_xp), dir_yp) * dxp * dyp; 
	rho_Z_av = sum(sum(rho, 1), 2) * dX * dY; 
	rho_Z_av = rho_Z_av(:);
	Z = (0 : NZ-1) * dZ;
	Z_bc = [-S.FDn:-1,NZ:NZ+S.FDn-1] * dZ;

	% a matrix of size Nz' x Nz', a(i,j) =  z_i - z_j
	%ZZp = abs(Z_bc.' - Z); 
	ZZp = abs(colminusrow(Z_bc.', Z)); 

	% V_av(Z) = int (rho_Z_av * |Z - Z'|) dZ', note this can be simplified
	% for charge neutrual systems to int (rho_Z_av * Z') dZ', indep of Z
	V_Z_av = (-2*pi/A_XY*dZ) * sum(bsxfun(@times, rho_Z_av', ZZp),2);

	% wave vectors
	GX = ifftshift( -floor(NX/2) : ceil(NX/2)-1 ) * (2*pi/NX);
	GY = ifftshift( -floor(NY/2) : ceil(NY/2)-1 ) * (2*pi/NY);
	[GX3D, GY3D, Z3D] = ndgrid(GX, GY, Z);
	GR3D = sqrt(GX3D.^2 + GY3D.^2);
	V_hat = zeros(NX,NY,2*S.FDn);
	for k = 1:length(Z_bc)  % in fact, no need to calculate for all Z values
		V_hat(:,:,k) = V_hat(:,:,k) + (2*pi*dZ) * ...
			sum(rho_hat ./ GR3D .* exp(-abs(Z3D - Z_bc(k)) .* GR3D) , 3);
	end
	V_hat(1,1,:) = 0; % remove the 0 frequency term (constant term in real-space)
	
	II = (1+S.FDn):(NX+S.FDn);
	JJ = (1+S.FDn):(NY+S.FDn);

	% put phi back to the 3D matrix
	ind_z = [1:S.FDn,NZ+S.FDn+1:NZ+2*S.FDn];
	phi(II,JJ,ind_z) = reshape(kron(V_Z_av', ones(NX,NY)),NX,NY,length(Z_bc));

	% remove 0 before ifft2 to include more terms
	phi(II,JJ,ind_z) = phi(II,JJ,ind_z) + 0*ifft2(V_hat); 
	% phi(S.isIn) = 0; % any value inside the fd-grid will be set to 0

	% permute phi back to original directions
	phi = real(phi);
	phi = permute(phi,reverse_order);
elseif S.BC == 4
	% find boundary conditions for periodic in 1D, dirichlet in 2D
	% Calculate phi
	phi = zeros(S.Nx+2*S.FDn, S.Ny+2*S.FDn, S.Nz+2*S.FDn);

	cellsize  = [S.L1, S.L2, S.L3];
	gridsizes = [S.Nx, S.Ny, S.Nz];
	meshsizes = [S.dx, S.dy, S.dz];
	bcs       = [S.BCx,S.BCy,S.BCz];

	% find which direction has Periodic BC
	dir_Z = find(bcs == 0); 
	dir_X = mod(dir_Z, 3) + 1;
	dir_Y = mod(dir_X, 3) + 1;

	% once we find the direction, we assume that direction is the Z
	% direction, the other two directions are then called X, Y
	NX = gridsizes(dir_X); NY = gridsizes(dir_Y); NZ = gridsizes(dir_Z);
	LX = cellsize(dir_X);  LY = cellsize(dir_Y);  LZ = cellsize(dir_Z);
	dX = meshsizes(dir_X); dY = meshsizes(dir_Y); dZ = meshsizes(dir_Z);
	% A_XY = LX * LY; % area of (x',y') surface

	% reshape rho to 3D 
	rho = reshape(rho+S.b, S.Nx, S.Ny, S.Nz); % note here after rho = rho + b

	% permute rho so that the new Z direction has Periodic BC
	new_order = [dir_X, dir_Y, dir_Z]; % a permutation of [1,2,3]
	[~, reverse_order] = sort(new_order); % reverse order to get back
	rho = permute(rho,new_order);
	phi = permute(phi,new_order);
	% rho = permute(rho,reverse_order); % this will recover original

	% sum over Z direction, \int (\rho) dz / LZ
	rho_XY_av = sum(rho,3) * (dZ/LZ); % NX * NY
	rho_XY_av = rho_XY_av(:);
	X = (0 : NX-1) * dX;
	Y = (0 : NY-1) * dY;

	[XX,YY] = ndgrid(X,Y);

	I_ex = 1:NX+2*S.FDn;
	J_ex = 1:NY+2*S.FDn;
	[II_ex,JJ_ex] = ndgrid(I_ex,J_ex);

	% flag for points inside the domain
	% isIn = ones(NX+2*S.FDn,NY+2*S.FDn);
	% isIn(1:S.FDn,:) = 0;
	% isIn(:,1:S.FDn) = 0;
	% isIn(S.FDn+NX+1:end,:) = 0;
	% isIn(:,S.FDn+NY+1:end) = 0;
	isIn = zeros(NX+2*S.FDn,NY+2*S.FDn);
	isIn(S.FDn+1:NX+S.FDn, S.FDn+1:NY+S.FDn) = 1;

	% positions of nodes outside the domain
	isbc = find(~isIn);

	XX_bc = (II_ex(isbc)-1-S.FDn) * dX;
	YY_bc = (JJ_ex(isbc)-1-S.FDn) * dY;

	% ln((X-X')^2 + (Y-Y')^2)
	% ln_RRp = log((XX_bc - XX(:)').^2 + (YY_bc - YY(:)').^2 );
	ln_RRp = log(colminusrow(XX_bc,XX(:)').^2 + colminusrow(YY_bc,YY(:)').^2);
	
	% V_XY_av = int (-ln((X-X')^2 + (Y-Y')^2) * rho_XY_av) dX'dY'
	V_XY_av = (-dX*dY) * sum(bsxfun(@times, rho_XY_av', ln_RRp),2);

	V_XY_av_full = zeros(NX+2*S.FDn,NY+2*S.FDn);
	V_XY_av_full(isbc) = V_XY_av;

	phi = bsxfun(@plus, phi, V_XY_av_full);

	% permute phi back to original directions
	phi = real(phi);
	phi = permute(phi,reverse_order);
end

% Calculate boundary conditions, this part can be combined for 2D
% periodic, 1D periodic and all Dirichlet BC
dx2 = S.dx * S.dx;
dy2 = S.dy * S.dy;
dz2 = S.dz * S.dz;
d = zeros(S.Nx,S.Ny,S.Nz);

II = (1+S.FDn):(S.Nx+S.FDn);
JJ = (1+S.FDn):(S.Ny+S.FDn);
KK = (1+S.FDn):(S.Nz+S.FDn);

% only add charge correction on Dirichlet boundaries
if S.BCx == 1
	for p = 1:S.FDn
		d = d - S.w2(p+1)/dx2 * (phi(II+p,JJ,KK) + phi(II-p,JJ,KK));
	end
end
if S.BCy == 1
	for p = 1:S.FDn
		d = d - S.w2(p+1)/dy2 * (phi(II,JJ+p,KK) + phi(II,JJ-p,KK));
	end
end
if S.BCz == 1
	for p = 1:S.FDn
		d = d - S.w2(p+1)/dz2 * (phi(II,JJ,KK+p) + phi(II,JJ,KK-p));
	end
end
d = d(:);
f = f + d;

end



% tool function colminusrow
function xmy = colminusrow(x,y)
% A column vector x minus a row vector.
% In Matlab versions after R2018b, it's just x - y
if (size(x,2) ~= 1)
	error('ERROR: the first vector must be a column vector');
end
if (size(y,1) ~= 1)
	error('ERROR: the second vector must be a row vector');
end

[xx,yy] = ndgrid(x,y);
xmy = xx - yy;

end

