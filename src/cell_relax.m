function S = cell_relax(S)
% @brief    This function performs cell relaxation of with fixed fractional 
%           coordinates of the atoms. Currently only volume relaxation is
%           implemented.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech

% max dialation of volume allowed
% S.max_dilatation = 1.2;

% relax cell in the periodic dims
S.cellrelax_dims = [1-S.BCx, 1-S.BCy, 1-S.BCz];
S.cellrelax_ndim = sum(S.cellrelax_dims);

if S.cellrelax_ndim < 1
	error('Cell relaxation is not allowed for isolated systems!');
end

V = S.L1 * S.L2 * S.L3; % volume

lower_bound =  V / S.max_dilatation;
upper_bound =  V * S.max_dilatation;

[optVol,S,~] = Brent(@volume_relax_fun, S, lower_bound, upper_bound, 1e-4, S.TOL_RELAX_CELL, S.max_relax_it);

% cell size related to optimal volume
V_old = S.L1 * S.L2 * S.L3;
scal = nthroot(optVol / V_old, S.cellrelax_ndim);
if S.cellrelax_dims(1) == 1
	S.L1 = S.L1 * scal;
end
if S.cellrelax_dims(2) == 1
	S.L2 = S.L2 * scal;
end
if S.cellrelax_dims(3) == 1
	S.L3 = S.L3 * scal;
end
S.optVol = optVol; % save optimized volume

fprintf('\n');
fprintf(' *****************************************************************\n');
fprintf(' *  Optimized volume: %.9f  \n', optVol);
fprintf(' *  Optimized cell size: %.9f %.9f %.9f\n',S.L1,S.L2,S.L3);
fprintf(' *****************************************************************\n');
fprintf('\n');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [max_P_stress,S] = volume_relax_fun(S, Vol) 
% Calculate principle stress for new volume

% update cell length (mesh size) related initialization and atom positions
% TODO: Need to go through msparc.m and initialization.m to find out full
% list


S = reinitialize_cell_mesh(S, Vol);

fprintf('\n');
fprintf(' **************************************************\n');
fprintf(' Reinitialization for relax step # %d\n',S.Relax_iter);
fprintf(' Volume: %-15.10f\n',Vol);
fprintf(' Cell  : %.10f %.10f %.10f\n',S.L1,S.L2,S.L3);
fprintf(' Mesh  : %f %f %f\n',S.dx,S.dy,S.dz);
fprintf(' rb    : %.2f %.2f %.2f\n',[S.Atm(:).rb_x,S.Atm(:).rb_y,S.Atm(:).rb_z]');
fprintf(' **************************************************\n');

% evaluate stress
S.Calc_stress = 1; % make sure stress calculation is on
[~,~,S] = electronicGroundStateAtomicForce(reshape(S.Atoms',[],1),S);

% print the result in .geopt file
if S.PrintRelaxout == 1
	PrintCellRelax(S);
end

S.ForceCount = S.ForceCount + 1;
S.Relax_iter = S.Relax_iter + 1;

% find the max principle stress in the relaxed dimensions
ind = (S.cellrelax_dims ~= 0);

%max_P_stress = max(eig(S.Stress(ind,ind)));
max_P_stress = sum(diag(S.Stress(ind,ind)))/S.cellrelax_ndim;

end


function S = reinitialize_cell_mesh(S, Vol)
% Once Volume is updated, the cell dimensions and mesh sized will be
% updated. All variables related to cell length and mesh size have to be
% updated as well.
%
% S.L1, S.L2, S.L3
% S.dx, S.dy, S.dz, S.dV
% finite-difference weights
% S.Atoms
% S.rb
% discrete Laplacian, gradient

if Vol <= 0.0
	error('Volume has become 0 or negative!');
end

Vol_old = S.L1 * S.L2 * S.L3; % old volume

% scaling factor
scal = nthroot(Vol / Vol_old, S.cellrelax_ndim);

% Calculate new length corresponding to volume Vol
if S.cellrelax_dims(1) == 1
	S.L1 = S.L1 * scal;
end
if S.cellrelax_dims(2) == 1
	S.L2 = S.L2 * scal;
end
if S.cellrelax_dims(3) == 1
	S.L3 = S.L3 * scal;
end

% Calculate new atom Cartesian coordinates
if S.cellrelax_dims(1) == 1
	S.Atoms(:,1) = S.Atoms(:,1) * scal;
	for ityp = 1:S.n_typ
		if (S.IsFrac(ityp) == 0)
			S.Atm(ityp).coords(:,1) = S.Atm(ityp).coords(:,1) * scal;
		end
	end
end

if S.cellrelax_dims(2) == 1
	S.Atoms(:,2) = S.Atoms(:,2) * scal;
	for ityp = 1:S.n_typ
		if (S.IsFrac(ityp) == 0)
			S.Atm(ityp).coords(:,2) = S.Atm(ityp).coords(:,2) * scal;
		end
	end
end

if S.cellrelax_dims(3) == 1
	S.Atoms(:,3) = S.Atoms(:,3) * scal;
	for ityp = 1:S.n_typ
		if (S.IsFrac(ityp) == 0)
			S.Atm(ityp).coords(:,3) = S.Atm(ityp).coords(:,3) * scal;
		end
	end
end


% update mesh size
if S.cellrelax_dims(1) == 1
	S.dx = S.dx * scal;
end
if S.cellrelax_dims(2) == 1
	S.dy = S.dy * scal;
end
if S.cellrelax_dims(3) == 1
	S.dz = S.dz * scal;
end
S.dV = S.dx * S.dy * S.dz * S.Jacb;

% Weights for spatial integration over domain
% Weights for spatial integration over domain
% S.W = IntgWts(S.Nx,S.Ny,S.Nz,S.BCx,S.BCy,S.BCz,S.xin,S);
if S.cell_typ == 1 || S.cell_typ == 2
	S.W = ones(S.N,1) * (S.dx*S.dy*S.dz*S.Jacb);
else
	S.W = IntgWts(S.Nx,S.Ny,S.Nz,S.BCx,S.BCy,S.BCz,S.xin,S);
end
% Brillouin-Zone Sampling
S = Generate_kpts(S);

% Create spherical harmonics for poisson solve for isolated clusters
if (S.BCx == 1 && S.BCy == 1 && S.BCz == 1)
	% Calculate Spherical Harmonics with origin shifted to the center of the domain
	xx_aug = (0-S.FDn:S.Nx+S.FDn-1)*S.dx;% - L1/2;
	yy_aug = (0-S.FDn:S.Ny+S.FDn-1)*S.dy;% - L2/2;
	zz_aug = (0-S.FDn:S.Nz+S.FDn-1)*S.dz;% - L3/2;
	[XX_AUG_3D,YY_AUG_3D,ZZ_AUG_3D] = ndgrid(xx_aug,yy_aug,zz_aug);
	% Find distances
	RR_AUG_3D = calculateDistance(XX_AUG_3D,YY_AUG_3D,ZZ_AUG_3D,S.L1/2,S.L2/2,S.L3/2,S);

	XX_AUG = reshape(XX_AUG_3D,[],1);
	YY_AUG = reshape(YY_AUG_3D,[],1);
	ZZ_AUG = reshape(ZZ_AUG_3D,[],1);
	RR_AUG = reshape(RR_AUG_3D,[],1);

	S.RR_AUG = RR_AUG;
	S.RR_AUG_3D = RR_AUG_3D;

	in_flag = ones(S.Nx+2*S.FDn,S.Ny+2*S.FDn,S.Nz+2*S.FDn);
	in_flag(1:S.FDn,:,:) = 0;
	in_flag(:,1:S.FDn,:) = 0;
	in_flag(:,:,1:S.FDn) = 0;
	in_flag(S.FDn+S.Nx+1:end,:,:) = 0;
	in_flag(:,S.FDn+S.Ny+1:end,:) = 0;
	in_flag(:,:,S.FDn+S.Nz+1:end) = 0;
	isIn = (in_flag~=0);

	S.isIn = isIn;
	S.RR = RR_AUG(isIn);

	pos_node_cart = coordinateTransformation(S,[XX_AUG,YY_AUG,ZZ_AUG],'noncart2cart_dis');
	pos_atm_cart = coordinateTransformation(S,[S.L1/2,S.L2/2,S.L3/2],'noncart2cart_dis');
	XX_AUG = bsxfun(@minus,pos_node_cart(:,1),pos_atm_cart(:,1));
	YY_AUG = bsxfun(@minus,pos_node_cart(:,2),pos_atm_cart(:,2));
	ZZ_AUG = bsxfun(@minus,pos_node_cart(:,3),pos_atm_cart(:,3));

	l_cut = 6;
	SH = repmat(struct([]),l_cut+1,1);
	for l = 0:l_cut
		for m = -l:l
			SH(l+1).Ylm_AUG(:,m+l+1) = sphericalHarmonics(XX_AUG,YY_AUG,ZZ_AUG,l,m,'real');
			Ylm_AUG_TEMP = SH(l+1).Ylm_AUG(:,m+l+1);
			SH(l+1).Ylm(:,m+l+1) = Ylm_AUG_TEMP(isIn);
			% SH(l+1).Ylm(:,m+l+1) = sphericalHarmonics(XX,YY,ZZ,l,m,'real');
		end
	end

	S.l_cut = l_cut;
	S.SH = SH;
end

% first find effective mesh size
dx2_inv = 1/(S.dx * S.dx);
dy2_inv = 1/(S.dy * S.dy);
dz2_inv = 1/(S.dz * S.dz);
h_eff = sqrt(3.0 / (dx2_inv + dy2_inv + dz2_inv));

% Update chebyshev polynomial degree if larger npl is required
S.npl =  max(Mesh2ChebDegree(h_eff), S.npl); 

% update precond_tol if more strict tol is required for new mesh size
S.precond_tol = min(h_eff * h_eff * 0.001, S.precond_tol);

% Calculate rb
S = Calculate_rb(S);

% Calculate discrete laplacian (1D) and discrete gradient indices' values
S = lapIndicesValues_1d(S);
S = gradIndicesValues(S);

% Calculate discrete laplacian
[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,[0 0 0]);
if S.cell_typ < 3
	S.Lap_std = S.lapc_T(1,1) * kron(speye(S.Nz),kron(speye(S.Ny),DL11))  +  S.lapc_T(2,2) * kron(speye(S.Nz),kron(DL22,speye(S.Nx))) + ...
				S.lapc_T(3,3) * kron(DL33,kron(speye(S.Ny),speye(S.Nx))) ;
	if (S.cell_typ == 2)
		MDL = S.lapc_T(1,2) * kron(speye(S.Nz),kron(DG2,DG1))  +  S.lapc_T(2,3) * kron(DG3,kron(DG2,speye(S.Nx))) + ...
			  S.lapc_T(1,3) * kron(DG3,kron(speye(S.Ny),DG1)) ;
		S.Lap_std = S.Lap_std + MDL;
	end
elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
	S.Lap_std = kron(speye(S.Nz),kron(speye(S.Ny),(DL11+DG1))) + kron(DL33,kron(speye(S.Ny),speye(S.Nx))) + ...
				kron(speye(S.Nz),kron(DL22,S.R2inv));
	if (S.cell_typ == 4 || S.cell_typ == 5)
		MDL = kron(speye(S.Nz),kron(DL22,speye(S.Nx))) +   kron(DG3,kron(DG2,speye(S.Nx)));
		S.Lap_std = S.Lap_std + MDL;
	end
end

% Calculate discrete gradient
S.grad_1 = blochGradient(S,[0 0 0],1);
S.grad_2 = blochGradient(S,[0 0 0],2);	
S.grad_3 = blochGradient(S,[0 0 0],3);

% Calculate preconditioners for negative discrete laplacian
[S.LapPreconL, S.LapPreconU] = ilu(S.Lap_std,struct('droptol',1e-5));

end


function [x, S, err] = Brent(FUN, S, x1, x2, tol, tol_f, max_iter) 
%@brief   Brent's method for finding the root of a function.
%
%@param FUN(x;S)  Function handle who's root will be evaluated.
%@param tol       Tolerance for minimum x interval to perform bisection.
%@param tol_f     Tolerance for f(x) to be considered close enough to 0
%
%@ref     W.H. Press, Numerical recepies 3rd edition: The art of scientific 
%          computing, Cambridge university press, 2007.
%
EPSILON = eps; % double precision 2^-52 = 2.22e-16
err = 0;

a = x1;
b = x2;
c = x2;

[fa,S] = FUN(S,a); 
[fb,S] = FUN(S,b);

%printf("f(%f) = %f, f(%f) = %f\n", x1, fa, x2, fb);
if((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))   
	x = 0.0; err = 1; % error flag
	error('Root must be within the given range (x1, x2) in Brent''s method');
end
fc = fb;

for iter = 1:max_iter
	if((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
		c = a; fc = fa;
		e = b-a;
		d = b-a;
	end
	if ( abs(fc) < abs(fb) )
		a = b; b = c; c = a;
		fa = fb; fb = fc; fc = fa;
	end
	tol1 = 2.0 * EPSILON * abs(b) + 0.5 * tol; 
	xm = 0.5*(c-b);
	%if(abs(xm) <= tol1 || abs(fb) < EPSILON)
	if(abs(xm) <= tol1 || abs(fb) < tol_f)
		x = b;
		fprintf(' Brent''s (bisection) method converged in %d iterations\n',iter);
		return;
	end
	if(abs(e) >= tol1 && abs(fa) > abs(fb))
		% attempt inverse quadratic interpolation
		s = fb / fa;
		if( a == c)
			p = 2.0 * xm * s;
			q = 1.0 - s;
		else
			q = fa / fc; r = fb / fc;
			p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
			q = (q - 1.0) * (r - 1.0) * (s - 1.0);
		end
		if(p > 0.0)
			%check whether in bounds
			q = -q;
		end
		p = abs(p);
		tol1q = tol1 * q;
		min1 = 3.0 * xm * q - abs(tol1q);
		eq = e * q;
		min2 = abs(eq);
		if(2.0 * p < min(min1,min2)) 
			%accept interpolation
			e = d; d = p / q;
		else
			% Bounds decreasing too slowly, use bisection
			d = xm; e = d;
		end
	else
		d = xm; e = d;
	end
	% move last best guess to a
	a = b; fa = fb;

	if (abs(d) > tol1)
		% evaluate new trial root
		b = b + d;
	else
		b = b + SIGN(tol1, xm);	
	end
	[fb,S] = FUN(S,b);
end

%x = 0.0; 
err = 1; % error flag
fprintf(' Maximum iterations exceeded in brents root finding method...exiting\n');

end


%%%%%%%%%%%%%%%%%%%%%%%
function y = SIGN(a,b)    
% y = abs(a) * sign(b);
if b > 0.0
	y = abs(a);
else
	y = -abs(a);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function npl = Mesh2ChebDegree(h) 
	% the relation between h and npl is fit with a cubic polynomial
	% p(x) = p3 * x^3 + p2 * x^2 + p1 * x + p0.
	p3 = -700. / 3.;
	p2 = 1240. / 3.;
	p1 = -773. / 3.;
	p0 = 1078. / 15.;
	if (h > 0.7) 
		npl = 14;
	else 
		npl = ((p3 * h + p2) * h + p1) * h + p0;
	end
	npl = round(npl);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = Calculate_rb(S)
% Starting and ending indices of b-region
if (S.cell_typ == 1 || S.cell_typ == 2)
	pos_atm_x = 0; % atom location in x-direction
	pos_atm_y = 0; % atom location in y-direction
	pos_atm_z = 0; % atom location in z-direction
	rb_up_x = 12;
	f_rby = @(y) y;
	rb_up_y = f_rby(12);
	rb_up_z = 12;
	
elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
	pos_atm_x = S.xmax_at; % maximum R coordinate of any atom
	pos_atm_y = 0; % atom location in theta-direction
	pos_atm_z = 0; % atom location in z-direction
	rb_up_x = S.xvac; % Radial direction vacuum
	f_rby = @(y) acos(1 - y^2/(2*pos_atm_x^2));
	rb_up_y = f_rby(12); % Theta direction
	rb_up_z = 12; % z-direction
	
end
ii_s_temp = -ceil(rb_up_x/S.dx);
ii_e_temp = ceil(rb_up_x/S.dx);
jj_s_temp = -ceil(rb_up_y/S.dy);
jj_e_temp = ceil(rb_up_y/S.dy);
kk_s_temp = 0;
kk_e_temp = ceil(rb_up_z/S.dz);
xx_temp = pos_atm_x + (ii_s_temp-S.FDn:ii_e_temp+S.FDn)*S.dx;
yy_temp = pos_atm_y + (jj_s_temp-S.FDn:jj_e_temp+S.FDn)*S.dy;
zz_temp = pos_atm_z + (kk_s_temp-S.FDn:kk_e_temp+S.FDn)*S.dz;
[XX_3D_temp,YY_3D_temp,ZZ_3D_temp] = ndgrid(xx_temp,yy_temp,zz_temp);
Nx = (ii_e_temp-ii_s_temp)+1;
Ny = (jj_e_temp-jj_s_temp)+1;
Nz = (kk_e_temp-kk_s_temp)+1;
% Find distances
dd_temp = calculateDistance(XX_3D_temp,YY_3D_temp,ZZ_3D_temp,pos_atm_x,pos_atm_y,pos_atm_z,S);

% Find integration weights
W_temp = IntgWts(Nx,Ny,Nz,1,1,1,xx_temp(S.FDn+1),S); % 1 - dirichlet BC on the boundary nodes
W_temp = reshape(W_temp,Nx,Ny,Nz);

% Find VJ and bJ
for ityp = 1:S.n_typ
	V_PS_temp = zeros(size(dd_temp));
	IsLargeThanRmax = dd_temp > S.Atm(ityp).r_grid_vloc(end);
	V_PS_temp(IsLargeThanRmax) = -S.Atm(ityp).Z;
	V_PS_temp(~IsLargeThanRmax) = interp1(S.Atm(ityp).r_grid_vloc, ...
		S.Atm(ityp).r_grid_vloc.*S.Atm(ityp).Vloc, dd_temp(~IsLargeThanRmax), 'spline');
	
	V_PS_temp = V_PS_temp./dd_temp;
	V_PS_temp(dd_temp<S.Atm(ityp).r_grid_vloc(2)) = S.Atm(ityp).Vloc(1);
	II_temp = 1+S.FDn : size(V_PS_temp,1)-S.FDn;
	JJ_temp = 1+S.FDn : size(V_PS_temp,2)-S.FDn;
	KK_temp = 1+S.FDn : size(V_PS_temp,3)-S.FDn;
	
	b_temp = pseudochargeDensity_atom(V_PS_temp,II_temp,JJ_temp,KK_temp,xx_temp(1),S);
	b_temp = -b_temp / (4*pi);
	err_rb = 100;
	count = 1;
	rb_x = S.Atm(ityp).rc;
	rb_y = f_rby(S.Atm(ityp).rc);
	rb_z = S.Atm(ityp).rc;
    rb_x = ceil(rb_x/S.dx-1e-12)*S.dx;
    rb_y = ceil(rb_y/S.dy-1e-12)*S.dy;
    rb_z = ceil(rb_z/S.dz-1e-12)*S.dz;
	fprintf(' Finding rb for %s ...\n',S.Atm(ityp).typ);
	while (err_rb > S.pseudocharge_tol && count <= 100 && rb_x <= rb_up_x && rb_y <= rb_up_y && rb_z <= rb_up_z )
		rb_x = rb_x + S.dx;
		rb_z = rb_z + S.dz;
		rb_y = f_rby(max(rb_x,rb_z));
		ii_rb = -1*ii_s_temp+S.FDn-floor(rb_x/S.dx)+1:-1*ii_s_temp+S.FDn+floor(rb_x/S.dx)+1;
		jj_rb = -1*jj_s_temp+S.FDn-floor(rb_y/S.dy)+1:-1*jj_s_temp+S.FDn+floor(rb_y/S.dy)+1;
		kk_rb = S.FDn+1:S.FDn+floor(rb_z/S.dz)+1; 
		err_rb = abs(sum(sum(sum(W_temp(ii_rb-S.FDn,jj_rb-S.FDn,kk_rb-S.FDn).*b_temp(ii_rb,jj_rb,kk_rb))))*2 + S.Atm(ityp).Z);
		fprintf(' rb = {%.3f %.3f %.3f}, int_b = %.15f, err_rb = %.3e\n',rb_x,rb_y,rb_z,2*sum(sum(sum(W_temp(ii_rb-S.FDn,jj_rb-S.FDn,kk_rb-S.FDn).*b_temp(ii_rb,jj_rb,kk_rb)))),err_rb);
		count = count + 1;
	end
	
	assert(rb_x<=rb_up_x && rb_y<=rb_up_y && rb_z<=rb_up_z,'Need to increase upper bound for rb!');
	S.Atm(ityp).rb_x = rb_x;
	S.Atm(ityp).rb_y = rb_y;
	S.Atm(ityp).rb_z = rb_z;
	% S.Atm(ityp).rb_x = ceil(rb_x/S.dx-1e-12)*S.dx; % + S.dx;
	% S.Atm(ityp).rb_y = ceil(rb_y/S.dy-1e-12)*S.dy; % + S.dy;
	% S.Atm(ityp).rb_z = ceil(rb_z/S.dz-1e-12)*S.dz; % + S.dz;
    fprintf(' rb = {%.3f %.3f %.3f}\n',S.Atm(ityp).rb_x,S.Atm(ityp).rb_y,S.Atm(ityp).rb_z);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W] = IntgWts(Nx,Ny,Nz,BCx,BCy,BCz,xin,S)
	if S.cell_typ == 1 || S.cell_typ == 2
		W_x = ones(Nx,1)*S.dx;
		W_x(1) = W_x(1) * (1-BCx*0.5);
		W_x(Nx) = W_x(Nx) * (1-BCx*0.5);

		W_y = ones(Ny,1)*S.dy;
		W_y(1) = W_y(1) * (1-BCy*0.5);
		W_y(Ny) = W_y(Ny) * (1-BCy*0.5);

		W_z = ones(Nz,1)*S.dz;
		W_z(1) = W_z(1) * (1-BCz*0.5);
		W_z(Nz) = W_z(Nz) * (1-BCz*0.5);

		W = kron(W_z,kron(W_y,W_x)) * S.Jacb;
	elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
		W = IntgWts_cychel(Nx,Ny,Nz,BCx,BCy,BCz,xin,S);
	end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrintCellRelax(S)
% @ brief   Function to write relaxation output in .cellopt file
%================================================================     
	if S.Relax_iter == 1
		fileID = fopen(S.relaxfname,'w');
	else
		fileID = fopen(S.relaxfname,'a');
	end
	fprintf(fileID,':RELAXSTEP: %d\n',S.Relax_iter);
	
	fprintf(fileID,':CELL: %18.10E %18.10E %18.10E\n',S.L1,S.L2,S.L3);
	fprintf(fileID,':VOLUME: %18.10E\n',S.L1*S.L2*S.L3);
	fprintf(fileID,':LATVEC:\n');
	fprintf(fileID,'%18.10E %18.10E %18.10E \n',S.lat_vec(1,:));
	fprintf(fileID,'%18.10E %18.10E %18.10E \n',S.lat_vec(2,:));
	fprintf(fileID,'%18.10E %18.10E %18.10E \n',S.lat_vec(3,:));
	fprintf(fileID,':STRESS:\n');
	fprintf(fileID,'%18.10E %18.10E %18.10E \n',S.Stress(1,:));
	fprintf(fileID,'%18.10E %18.10E %18.10E \n',S.Stress(2,:));
	fprintf(fileID,'%18.10E %18.10E %18.10E \n',S.Stress(3,:));
	
	fclose(fileID);
end




function [S] = Generate_kpts(S)
	nkpt = S.nkpt;
	if (S.BCx == 1 && nkpt(1) > 1)
		error(' nkpt cannot be greater than 1 in Dirichlet boundary direction (x)');
	end
	if (S.BCy == 1 && nkpt(2) > 1)
		error(' nkpt cannot be greater than 1 in Dirichlet boundary direction (y)');
	end
	if (S.BCz == 1 && nkpt(3) > 1)
		error(' nkpt cannot be greater than 1 in Dirichlet boundary direction (z)');
	end

	% Monkhorst-pack grid for Brillouin zone sampling
	MPG_typ1 = @(nkpt) (2*(1:nkpt) - nkpt - 1)/2; % MP grid points for infinite group order
	MPG_typ2 = @(nkpt) (0:nkpt-1); % MP grid points for finite group order

	if S.cell_typ < 3
		kptgrid_x = (1/nkpt(1)) * MPG_typ1(nkpt(1));
		kptgrid_y = (1/nkpt(2)) * MPG_typ1(nkpt(2));
		kptgrid_z = (1/nkpt(3)) * MPG_typ1(nkpt(3));
		sumx = 0;
		sumy = 0; 
		sumz = 0;
		% shift kpoint grid 
		kptgrid_x = kptgrid_x + S.kptshift(1) * (1/nkpt(1));
		kptgrid_y = kptgrid_y + S.kptshift(2) * (1/nkpt(2));
		kptgrid_z = kptgrid_z + S.kptshift(3) * (1/nkpt(3));
	
		% map k-points back to BZ
		temp_epsilon = eps; % include the right boundary k-points instead of left
		kptgrid_x = mod(kptgrid_x + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
		kptgrid_y = mod(kptgrid_y + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
		kptgrid_z = mod(kptgrid_z + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
	elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
		kptgrid_x = (1/nkpt(1)) * MPG_typ1(nkpt(1));
		kptgrid_y = (1/nkpt(2)) * MPG_typ2(nkpt(2));
		kptgrid_z = (1/nkpt(3)) * MPG_typ1(nkpt(3));
		sumx = 0;
		sumy = nkpt(2); 
		sumz = 0;
	end    
	
	% Scale kpoints
	kptgrid_x = (2*pi/S.L1) * kptgrid_x;
	kptgrid_y = (2*pi/S.L2) * kptgrid_y;
	kptgrid_z = (2*pi/S.L3) * kptgrid_z;

	[kptgrid_X, kptgrid_Y, kptgrid_Z] = ndgrid(kptgrid_x,kptgrid_y,kptgrid_z);
	kptgrid = [reshape(kptgrid_X,[],1),reshape(kptgrid_Y,[],1),reshape(kptgrid_Z,[],1)];
	disp('kpoint grid before symmetry:');
	disp(kptgrid);
	
	tnkpt = prod(nkpt);
	wkpt = ones(tnkpt,1)/tnkpt;% weights for k-points
	TOL = 1e-8;
	% Time-Reversal Symmetry to reduce k-points
	if S.TimeRevSym == 1
		Ikpt = zeros(tnkpt,1);
		Ikpt_rev = zeros(tnkpt,1);
		for ii = 1:tnkpt
			for jj = ii+1:tnkpt
				if (abs(kptgrid(ii,1) + kptgrid(jj,1) - sumx) < TOL) && (abs(kptgrid(ii,2) + kptgrid(jj,2) - sumy) < TOL) && (abs(kptgrid(ii,3) + kptgrid(jj,3) - sumz) < TOL)
					Ikpt(ii) = 1;
					Ikpt_rev(jj) = 1;
				end
			end
		end
		Ikpt = Ikpt>0.5;
		Ikpt_rev = Ikpt_rev>0.5;
		wkpt(Ikpt_rev) = 2*wkpt(Ikpt_rev);
		kptgrid = kptgrid(~Ikpt,:);
		wkpt = wkpt(~Ikpt);
		tnkpt = size(wkpt,1);
	end

	disp('kpoint grid after symmetry:');	
	disp(kptgrid);
	% Store into the structure
	S.kptgrid = kptgrid;
	S.tnkpt   = tnkpt;
	S.wkpt    = wkpt;
end


