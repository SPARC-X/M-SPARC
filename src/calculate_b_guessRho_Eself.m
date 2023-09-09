function S = calculate_b_guessRho_Eself(S)
% @brief    calculate_b_guessRho_Eself(S) calculates the pseudocharge (& ref), 
%           self energy (& ref), electronic density guess, electrostatic 
%           energy correction and electrostatic potential correction.
%
% @param S  A struct that contains the relevant fields.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%

t1 = tic;
fprintf('\n Starting pseudocharge generation and self energy calculation...\n');

% initialization
S.b = zeros(S.N,1);
S.b_ref = zeros(S.N,1);
S.rho_at = zeros(S.N,1);
if S.spin_typ == 1
	mz = zeros(S.N,1);
elseif S.spin_typ == 2
    mx = zeros(S.N,1);
    my = zeros(S.N,1);
    mz = zeros(S.N,1);
end
S.Eself = 0;
S.Eself_ref = 0;
S.V_c = zeros(S.N,1);

% Pseudocharge generation and self energy calculation
count_typ = 1;
count_typ_atms = 1;
for JJ_a = 1:S.n_atm % loop over all the atoms
	% Atom position of atom JJ_a
	x0 = S.Atoms(JJ_a,1);
	y0 = S.Atoms(JJ_a,2);
	z0 = S.Atoms(JJ_a,3);
	
	% Note the S.dx, S.dy, S.dz terms are to ensure the image rb-region overlap w/ fund. domain
	if S.BCx == 0
		n_image_xl = floor((S.Atoms(JJ_a,1) + S.Atm(count_typ).rb_x)/S.L1);
		n_image_xr = floor((S.L1 - S.Atoms(JJ_a,1)+S.Atm(count_typ).rb_x-S.dx)/S.L1);
	else
		n_image_xl = 0;
		n_image_xr = 0;
	end
	
	if S.BCy == 0
		n_image_yl = floor((S.Atoms(JJ_a,2) + S.Atm(count_typ).rb_y)/S.L2);
		n_image_yr = floor((S.L2 - S.Atoms(JJ_a,2)+S.Atm(count_typ).rb_y-S.dy)/S.L2);
	else
		n_image_yl = 0;
		n_image_yr = 0;
	end
	
	if S.BCz == 0
		n_image_zl = floor((S.Atoms(JJ_a,3) + S.Atm(count_typ).rb_z)/S.L3);
		n_image_zr = floor((S.L3 - S.Atoms(JJ_a,3)+S.Atm(count_typ).rb_z-S.dz)/S.L3);
	else
		n_image_zl = 0;
		n_image_zr = 0;
	end
	
	% Total No. of images of atom JJ_a (including atom JJ_a)
	n_image_total = (n_image_xl+n_image_xr+1) * (n_image_yl+n_image_yr+1) * (n_image_zl+n_image_zr+1);
	% Find the coordinates for all the images
	xx_img = (-n_image_xl : n_image_xr) * S.L1 + x0;
	yy_img = (-n_image_yl : n_image_yr) * S.L2 + y0;
	zz_img = (-n_image_zl : n_image_zr) * S.L3 + z0;
	[XX_IMG_3D,YY_IMG_3D,ZZ_IMG_3D] = ndgrid(xx_img,yy_img,zz_img);
	
	% Loop over all image(s) of atom JJ_a (including atom JJ_a)
	for count_image = 1:n_image_total
		% Atom position of the image
		x0_i = XX_IMG_3D(count_image);
		y0_i = YY_IMG_3D(count_image);
		z0_i = ZZ_IMG_3D(count_image);
		
		%************************************************************************
		%*          Calculate b, b_ref, Eself, Eself_ref and rho_at             *
		%************************************************************************
		% Starting and ending indices of b-region
		ii_s = ceil ((x0_i - S.Atm(count_typ).rb_x)/S.dx) + 1;
		ii_e = floor((x0_i + S.Atm(count_typ).rb_x)/S.dx) + 1;
		jj_s = ceil ((y0_i - S.Atm(count_typ).rb_y)/S.dy) + 1;
		jj_e = floor((y0_i + S.Atm(count_typ).rb_y)/S.dy) + 1;
		kk_s = ceil ((z0_i - S.Atm(count_typ).rb_z)/S.dz) + 1;
		kk_e = floor((z0_i + S.Atm(count_typ).rb_z)/S.dz) + 1;
		
		% Check if the b-region is inside the domain in Dirichlet BC
		% direction
		isInside = (S.BCx == 0 || (S.BCx == 1 && (ii_s>1) && (ii_e<S.Nx))) && ...
		   (S.BCy == 0 || (S.BCy == 1 && (jj_s>1) && (jj_e<S.Ny))) && ...
		   (S.BCz == 0 || (S.BCz == 1 && (kk_s>1) && (kk_e<S.Nz)));
		% assert(isInside,'ERROR: Atom too close to boundary for b calculation');
		if ~isInside
			fprintf(' WARNING: Atom %d too close to boundary for b calculation\n',JJ_a);
		end
		ii_s = max(ii_s,1);
		ii_e = min(ii_e,S.Nx);
		jj_s = max(jj_s,1);
		jj_e = min(jj_e,S.Ny);
		kk_s = max(kk_s,1);
		kk_e = min(kk_e,S.Nz);

		xx = S.xin + (ii_s-S.FDn-1:ii_e+S.FDn-1)*S.dx;
		yy = S.yin + (jj_s-S.FDn-1:jj_e+S.FDn-1)*S.dy;
		zz = S.zin + (kk_s-S.FDn-1:kk_e+S.FDn-1)*S.dz;
		[XX_3D,YY_3D,ZZ_3D] = ndgrid(xx,yy,zz);

		% Find distances
		dd = calculateDistance(XX_3D,YY_3D,ZZ_3D,x0_i,y0_i,z0_i,S);
		
		% Pseudopotential at grid points through interpolation
		V_PS = zeros(size(dd));
		IsLargeThanRmax = dd > S.Atm(count_typ).r_grid_vloc(end);
		V_PS(IsLargeThanRmax) = -S.Atm(count_typ).Z;
		V_PS(~IsLargeThanRmax) = interp1(S.Atm(count_typ).r_grid_vloc, S.Atm(count_typ).r_grid_vloc.*S.Atm(count_typ).Vloc, dd(~IsLargeThanRmax), 'spline');

		V_PS = V_PS./dd;
		V_PS(dd<S.Atm(count_typ).r_grid_vloc(2)) = S.Atm(count_typ).Vloc(1); % WARNING

		% Reference potential at grid points
		rc_ref = S.rc_ref;
		V_PS_ref = zeros(size(dd));
		I_ref = dd<rc_ref;
		V_PS_ref(~I_ref) = -(S.Atm(count_typ).Z)./dd(~I_ref);
		V_PS_ref(I_ref) = -S.Atm(count_typ).Z*(9*dd(I_ref).^7-30*rc_ref*dd(I_ref).^6 ...
			+28*rc_ref*rc_ref*dd(I_ref).^5-14*(rc_ref^5)*dd(I_ref).^2+12*rc_ref^7)/(5*rc_ref^8);
		
		% Isolated atom electron density at grid points through interpolation
		rho_isolated_atom = interp1(S.Atm(count_typ).r_grid_rho, S.Atm(count_typ).rho_isolated_guess, dd, 'spline');
		rho_isolated_atom(dd > S.Atm(count_typ).r_grid_rho(end)) = 0;
		
		% Pseudocharge density, rho_at, Eself calculation
		II = 1+S.FDn : size(V_PS,1)-S.FDn;
		JJ = 1+S.FDn : size(V_PS,2)-S.FDn;
		KK = 1+S.FDn : size(V_PS,3)-S.FDn;

		% Calculate bJ and bJ_ref
		bJ = pseudochargeDensity_atom(V_PS,II,JJ,KK,xx(1),S);
		bJ_ref = pseudochargeDensity_atom(V_PS_ref,II,JJ,KK,xx(1),S);
		
		[II_rb,JJ_rb,KK_rb] = ndgrid(ii_s:ii_e,jj_s:jj_e,kk_s:kk_e);
		Rowcount_rb = (KK_rb-1)*S.Nx*S.Ny + (JJ_rb-1)*S.Nx + II_rb;
		sz = size(S.b(Rowcount_rb));
		S.b(Rowcount_rb) = S.b(Rowcount_rb) + reshape(bJ(II,JJ,KK), sz);
		S.b_ref(Rowcount_rb) = S.b_ref(Rowcount_rb) + reshape(bJ_ref(II,JJ,KK),sz);
		rho_add = reshape(rho_isolated_atom(II,JJ,KK),sz);
		S.rho_at(Rowcount_rb) = S.rho_at(Rowcount_rb) + rho_add;
		if S.spin_typ == 1
            mz(Rowcount_rb) = mz(Rowcount_rb) + S.Atm(count_typ).mag(count_typ_atms,3)/S.Atm(count_typ).Z * rho_add;
        elseif S.spin_typ == 2
            mx(Rowcount_rb) = mx(Rowcount_rb) + S.Atm(count_typ).mag(count_typ_atms,1)/S.Atm(count_typ).Z * rho_add;
            my(Rowcount_rb) = my(Rowcount_rb) + S.Atm(count_typ).mag(count_typ_atms,2)/S.Atm(count_typ).Z * rho_add;
            mz(Rowcount_rb) = mz(Rowcount_rb) + S.Atm(count_typ).mag(count_typ_atms,3)/S.Atm(count_typ).Z * rho_add;
		end
		S.Eself = S.Eself + 0.5 * sum(sum(sum(bJ(II,JJ,KK).*V_PS(II,JJ,KK).*S.W(Rowcount_rb) )));
		S.Eself_ref = S.Eself_ref + 0.5 * sum(sum(sum(bJ_ref(II,JJ,KK).*V_PS_ref(II,JJ,KK).*S.W(Rowcount_rb) )));
		S.V_c(Rowcount_rb) = S.V_c(Rowcount_rb) + reshape(V_PS_ref(II,JJ,KK) - V_PS(II,JJ,KK),sz);
	end % end of loop over images (including atom JJ_a)
	
	% Check if same type of atoms are over
	if count_typ_atms == S.Atm(count_typ).n_atm_typ
		count_typ_atms = 1;
		count_typ = count_typ + 1;
	else
		count_typ_atms = count_typ_atms + 1;
	end
end % end of loop over atoms

% Scaling by factor -1/(4pi) for b
S.b = (-1/(4*pi))*S.b;
S.Eself = (-1/(4*pi))*S.Eself;

% Scaling by factor -1/(4pi) for b_ref
S.b_ref = (-1/(4*pi))*S.b_ref;
S.Eself_ref = (-1/(4*pi))*S.Eself_ref;

% find positive charges and negative charges
S.PosCharge = abs(dot(S.W, S.b));
S.NegCharge = -S.PosCharge + S.NetCharge;

% Check for accuracy of b
fprintf(' Integration b = %.12f\n\n',abs(dot(S.W,S.b)));
%assert(abs(dot(S.W,S.b)+S.Nelectron)/S.Nelectron < 1e-4,'Pseudocharge needs improved accuracy!')

% Check for accuracy of b_ref
fprintf(' Integration b_ref = %.12f\n\n',abs(dot(S.W,S.b_ref)));
%assert(abs(dot(S.W,S.b_ref)+S.Nelectron)/S.Nelectron < 1e-4,'Rerefence Pseudocharge needs improved accuracy!')


%*************************************************************
%*                    Designate rho_at                       *
%*************************************************************
% Guess for electron density
if S.spin_typ == 1
    S.mag = mz;
	S.rho_at = [S.rho_at 0.5*(S.rho_at+mz) 0.5*(S.rho_at-mz)];
elseif S.spin_typ == 2
    magnorm = sqrt(mx.^2 + my.^2 + mz.^2);
    S.mag = [magnorm mx my mz];
    S.rho_at = [S.rho_at 0.5*(S.rho_at+magnorm) 0.5*(S.rho_at-magnorm)];
end

rho_scal = abs(S.NegCharge/dot(S.W,S.rho_at(:,1)));
S.rho_at = rho_scal*S.rho_at;

if S.spin_typ ~= 0
	S.netM = sum(S.mag)*S.dV;
	fprintf('======================================\n');
    fprintf(' Net initial magnetization is: ');
    fprintf('%.6f ', S.netM);
    fprintf('\n');
	fprintf('======================================\n\n');
end

fprintf(' ****************************************\n');
fprintf(' *          Eself_ref = %f       *\n',S.Eself_ref);
fprintf(' ****************************************\n');

%*************************************************************
%*                     Calculate E_corr                      *
%*************************************************************
S.E_corr = 0.5*sum((S.b_ref+S.b).*S.V_c.*S.W) + S.Eself - S.Eself_ref;
%fprintf('Eself: %.15f, E_corr: %.15f \n',S.Eself,S.E_corr);

fprintf(' Done. (%f s)\n',toc(t1));


% Calculate core charge (NLCC) by summing all the atoms and their images for the fundamental domain
count_typ = 1;
count_typ_atms = 1;
S.rho_Tilde_at = zeros(S.Nx,S.Ny,S.Nz);
for JJ_a = 1:S.n_atm % loop over all the atoms
	% Atom position
	x0 = S.Atoms(JJ_a,1);
	y0 = S.Atoms(JJ_a,2);
	z0 = S.Atoms(JJ_a,3);
	% Note the S.dx, S.dy, S.dz terms are to ensure the image rb-region overlap w/ fund. domain
	if S.BCx == 0
		n_image_xl = floor((S.Atoms(JJ_a,1) + max(S.Atm(count_typ).rb_x))/S.L1);
		n_image_xr = floor((S.L1 - S.Atoms(JJ_a,1)+max(S.Atm(count_typ).rb_x))/S.L1);
	else
		n_image_xl = 0;
		n_image_xr = 0;
	end
	
	if S.BCy == 0
		n_image_yl = floor((S.Atoms(JJ_a,2) + max(S.Atm(count_typ).rb_y))/S.L2);
		n_image_yr = floor((S.L2 - S.Atoms(JJ_a,2)+max(S.Atm(count_typ).rb_y))/S.L2);
	else
		n_image_yl = 0;
		n_image_yr = 0;
	end
	
	if S.BCz == 0
		n_image_zl = floor((S.Atoms(JJ_a,3) + max(S.Atm(count_typ).rb_z))/S.L3);
		n_image_zr = floor((S.L3 - S.Atoms(JJ_a,3)+max(S.Atm(count_typ).rb_z))/S.L3);
	else
		n_image_zl = 0;
		n_image_zr = 0;
	end
	
	% Total No. of images of atom JJ_a (including atom JJ_a)
	n_image_total = (n_image_xl+n_image_xr+1) * (n_image_yl+n_image_yr+1) * (n_image_zl+n_image_zr+1);
	% Find the coordinates for all the images
	xx_img = [-n_image_xl : n_image_xr] * S.L1 + x0;
	yy_img = [-n_image_yl : n_image_yr] * S.L2 + y0;
	zz_img = [-n_image_zl : n_image_zr] * S.L3 + z0;
	[XX_IMG_3D,YY_IMG_3D,ZZ_IMG_3D] = ndgrid(xx_img,yy_img,zz_img);

	% Loop over all image(s) of atom JJ_a (including atom JJ_a)
	for count_image = 1:n_image_total

		% Atom position of the image
		x0_i = XX_IMG_3D(count_image);
		y0_i = YY_IMG_3D(count_image);
		z0_i = ZZ_IMG_3D(count_image);

		% Indices of closest grid point to atom
		pos_ii = round((x0_i-S.xin) / S.dx) + 1;
		pos_jj = round((y0_i-S.yin) / S.dy) + 1;
		pos_kk = round((z0_i-S.zin) / S.dz) + 1;

		% Starting and ending indices of b-region
		ii_s = pos_ii - ceil(max(S.Atm(count_typ).r_grid_rho_Tilde)/S.dx+0.5);
		ii_e = pos_ii + ceil(max(S.Atm(count_typ).r_grid_rho_Tilde)/S.dx+0.5);
		jj_s = pos_jj - ceil(max(S.Atm(count_typ).r_grid_rho_Tilde)/S.dy+0.5);
		jj_e = pos_jj + ceil(max(S.Atm(count_typ).r_grid_rho_Tilde)/S.dy+0.5);
		kk_s = pos_kk - ceil(max(S.Atm(count_typ).r_grid_rho_Tilde)/S.dz+0.5);
		kk_e = pos_kk + ceil(max(S.Atm(count_typ).r_grid_rho_Tilde)/S.dz+0.5);

		% Check if the b-region is inside the domain in Dirichlet BC
		% direction
		%isInside = (S.BCx == 0 || (S.BCx == 1 && (ii_s>1) && (ii_e<S.Nx))) && ...
		%   (S.BCy == 0 || (S.BCy == 1 && (jj_s>1) && (jj_e<S.Ny))) && ...
		%   (S.BCz == 0 || (S.BCz == 1 && (kk_s>1) && (kk_e<S.Nz)));
		% assert(isInside,'ERROR: Atom too close to boundary for b calculation');
		ii_s = max(ii_s,1);
		ii_e = min(ii_e,S.Nx);
		jj_s = max(jj_s,1);
		jj_e = min(jj_e,S.Ny);
		kk_s = max(kk_s,1);
		kk_e = min(kk_e,S.Nz);

		xx = S.xin + (ii_s-2*S.FDn-1:ii_e+2*S.FDn-1)*S.dx;% - x0_i;
		yy = S.yin + (jj_s-2*S.FDn-1:jj_e+2*S.FDn-1)*S.dy;% - y0_i;
		zz = S.zin + (kk_s-2*S.FDn-1:kk_e+2*S.FDn-1)*S.dz;% - z0_i;
    
		[XX_3D,YY_3D,ZZ_3D] = ndgrid(xx,yy,zz);

		% Find distances
		dd = calculateDistance(XX_3D,YY_3D,ZZ_3D,x0_i,y0_i,z0_i,S);

		% Pseudopotential at grid points through interpolation
		rho_Tilde_at= zeros(size(dd));
		IsLargeThanRmax = dd > max(S.Atm(count_typ).r_grid_rho_Tilde);
		rho_Tilde_at(IsLargeThanRmax) =0;
		rho_Tilde_at(~IsLargeThanRmax) = interp1(S.Atm(count_typ).r_grid_rho_Tilde,S.Atm(count_typ).rho_Tilde,  dd(~IsLargeThanRmax), 'spline');
		
		II = 1+2*S.FDn : size(rho_Tilde_at,1)-2*S.FDn;
		JJ = 1+2*S.FDn : size(rho_Tilde_at,2)-2*S.FDn;
		KK = 1+2*S.FDn : size(rho_Tilde_at,3)-2*S.FDn;        
		S.rho_Tilde_at(ii_s:ii_e,jj_s:jj_e,kk_s:kk_e) = S.rho_Tilde_at(ii_s:ii_e,jj_s:jj_e,kk_s:kk_e) + rho_Tilde_at(II,JJ,KK);
	end
	% Check if same type of atoms are over
	if count_typ_atms == S.Atm(count_typ).n_atm_typ
		count_typ_atms = 1;
		count_typ = count_typ + 1;
	else
		count_typ_atms = count_typ_atms + 1;
	end
end % end of loop over atoms
S.rho_Tilde_at = S.rho_Tilde_at(:);

end
