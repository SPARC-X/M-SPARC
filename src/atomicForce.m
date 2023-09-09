function force = atomicForce(S)
% @brief    atomicForce(S) calculates the atomic force.
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%==========================================================================

fprintf('\n Starting atomic force calculation ... \n');

% Calculate local forces
tic_locforces = tic;

% Dpseudo_x = S.grad_1*(S.b + S.b_ref);
% Dpseudo_y = S.grad_2*(S.b + S.b_ref);
% Dpseudo_z = S.grad_3*(S.b + S.b_ref);

Dphi_x = S.grad_1*(S.phi);
Dphi_y = S.grad_2*(S.phi);
Dphi_z = S.grad_3*(S.phi);

DVc_x = S.grad_1*S.V_c;
DVc_y = S.grad_2*S.V_c;
DVc_z = S.grad_3*S.V_c;

% initialization
force_corr = zeros(S.n_atm,3);
force_nloc = zeros(S.n_atm,3);
force_local = zeros(S.n_atm,3);
force_xc = zeros(S.n_atm,3);

if S.NLCC_flag
    count_typ = 1;
    count_typ_atms = 1;
    for JJ_a = 1:S.n_atm % loop over all the atoms
        % Atom position
        x0 = S.Atoms(JJ_a,1);
        y0 = S.Atoms(JJ_a,2);
        z0 = S.Atoms(JJ_a,3);
        % Note the S.dx, S.dy, S.dz terms are to ensure the image rb-region overlap w/ fund. domain
        if S.BCx == 0
            n_image_xl = floor((S.Atoms(JJ_a,1) + max(S.Atm(count_typ).r_grid_rho_Tilde))/S.L1);
            n_image_xr = floor((S.L1 - S.Atoms(JJ_a,1)+max(S.Atm(count_typ).r_grid_rho_Tilde))/S.L1);
        else
            n_image_xl = 0;
            n_image_xr = 0;
        end

        if S.BCy == 0
            n_image_yl = floor((S.Atoms(JJ_a,2) + max(S.Atm(count_typ).r_grid_rho_Tilde))/S.L2);
            n_image_yr = floor((S.L2 - S.Atoms(JJ_a,2)+max(S.Atm(count_typ).r_grid_rho_Tilde))/S.L2);
        else
            n_image_yl = 0;
            n_image_yr = 0;
        end

        if S.BCz == 0
            n_image_zl = floor((S.Atoms(JJ_a,3) + max(S.Atm(count_typ).r_grid_rho_Tilde))/S.L3);
            n_image_zr = floor((S.L3 - S.Atoms(JJ_a,3)+max(S.Atm(count_typ).r_grid_rho_Tilde))/S.L3);
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


            % calculating gradient of rho_Tilde
            drho_Tilde_at_x = zeros(size(rho_Tilde_at)); drho_Tilde_at_y = zeros(size(rho_Tilde_at)); drho_Tilde_at_z = zeros(size(rho_Tilde_at));

            II = 1+2*S.FDn : size(rho_Tilde_at,1)-2*S.FDn;
            JJ = 1+2*S.FDn : size(rho_Tilde_at,2)-2*S.FDn;
            KK = 1+2*S.FDn : size(rho_Tilde_at,3)-2*S.FDn;

            for p = 1:S.FDn
                drho_Tilde_at_x(II,JJ,KK) = drho_Tilde_at_x(II,JJ,KK) + S.w1(p+1)/S.dx*(rho_Tilde_at(II+p,JJ,KK)-rho_Tilde_at(II-p,JJ,KK));
                drho_Tilde_at_y(II,JJ,KK) = drho_Tilde_at_y(II,JJ,KK) + S.w1(p+1)/S.dy*(rho_Tilde_at(II,JJ+p,KK)-rho_Tilde_at(II,JJ-p,KK));
                drho_Tilde_at_z(II,JJ,KK) = drho_Tilde_at_z(II,JJ,KK) + S.w1(p+1)/S.dz*(rho_Tilde_at(II,JJ,KK+p)-rho_Tilde_at(II,JJ,KK-p));
            end
            % Calculate xc 2nd term force components
            [II_rb,JJ_rb,KK_rb] = ndgrid(ii_s:ii_e,jj_s:jj_e,kk_s:kk_e);
            Rowcount_rb = (KK_rb-1)*S.Nx*S.Ny + (JJ_rb-1)*S.Nx + II_rb;
            if S.spin_typ == 0
                force_xc(JJ_a,1) = force_xc(JJ_a,1) + sum(sum(sum( drho_Tilde_at_x(II,JJ,KK) .* ( S.Vxc(Rowcount_rb) ) .* S.W(Rowcount_rb) )));
                force_xc(JJ_a,2) = force_xc(JJ_a,2) + sum(sum(sum( drho_Tilde_at_y(II,JJ,KK) .* ( S.Vxc(Rowcount_rb) ) .* S.W(Rowcount_rb) )));
                force_xc(JJ_a,3) = force_xc(JJ_a,3) + sum(sum(sum( drho_Tilde_at_z(II,JJ,KK) .* ( S.Vxc(Rowcount_rb) ) .* S.W(Rowcount_rb) )));
            else
                vxc = S.Vxc(:,1)+S.Vxc(:,2);
                force_xc(JJ_a,1) = force_xc(JJ_a,1) + sum(sum(sum( 0.5*drho_Tilde_at_x(II,JJ,KK) .* (vxc(Rowcount_rb) ) .* S.W(Rowcount_rb) )));
                force_xc(JJ_a,2) = force_xc(JJ_a,2) + sum(sum(sum( 0.5*drho_Tilde_at_y(II,JJ,KK) .* (vxc(Rowcount_rb) ) .* S.W(Rowcount_rb) )));
                force_xc(JJ_a,3) = force_xc(JJ_a,3) + sum(sum(sum( 0.5*drho_Tilde_at_z(II,JJ,KK) .* (vxc(Rowcount_rb) ) .* S.W(Rowcount_rb) )));
            end
        end
        % Check if same type of atoms are over
        if count_typ_atms == S.Atm(count_typ).n_atm_typ
            count_typ_atms = 1;
            count_typ = count_typ + 1;
        else
            count_typ_atms = count_typ_atms + 1;
        end

    end % end of loop over atoms
end

count_typ = 1;
count_typ_atms = 1;
for JJ_a = 1:S.n_atm % loop over all the atoms   
	% Atom position
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
		
		%**********************************************************************
		%*                  Calculate local atomic force term                 *
		%**********************************************************************
		
		% Starting and ending indices of b-region
		ii_s = ceil ((x0_i - S.Atm(count_typ).rb_x)/S.dx) + 1;
		ii_e = floor((x0_i + S.Atm(count_typ).rb_x)/S.dx) + 1;
		jj_s = ceil ((y0_i - S.Atm(count_typ).rb_y)/S.dy) + 1;
		jj_e = floor((y0_i + S.Atm(count_typ).rb_y)/S.dy) + 1;
		kk_s = ceil ((z0_i - S.Atm(count_typ).rb_z)/S.dz) + 1;
		kk_e = floor((z0_i + S.Atm(count_typ).rb_z)/S.dz) + 1;
		
		 % Check if the b-region is inside the domain in Dirichlet BC
		 % direction
		isInside = (S.BCx == 0 || (S.BCx == 1 && (ii_s>=1) && (ii_e<=S.Nx))) && ...
			(S.BCy == 0 || (S.BCy == 1 && (jj_s>=1) && (jj_e<=S.Ny))) && ...
			(S.BCz == 0 || (S.BCz == 1 && (kk_s>=1) && (kk_e<=S.Nz)));
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
		V_PS(dd<S.Atm(count_typ).r_grid_vloc(2)) = S.Atm(count_typ).Vloc(1);
		
		% Reference potential at grid points
		rc_ref = S.rc_ref; % WARNING: Might need smaller if pseudocharges overlap
		V_PS_ref = zeros(size(dd));
		I_ref = dd<rc_ref;
		V_PS_ref(~I_ref) = -(S.Atm(count_typ).Z)./dd(~I_ref);
		V_PS_ref(I_ref) = -S.Atm(count_typ).Z*(9*dd(I_ref).^7-30*rc_ref*dd(I_ref).^6 ...
			+28*rc_ref*rc_ref*dd(I_ref).^5-14*(rc_ref^5)*dd(I_ref).^2+12*rc_ref^7)/(5*rc_ref^8);
		
		% Pseudocharge density
		II = 1+S.FDn : size(V_PS,1)-S.FDn;
		JJ = 1+S.FDn : size(V_PS,2)-S.FDn;
		KK = 1+S.FDn : size(V_PS,3)-S.FDn;
		
		% Calculate bJ and bJ_ref
		bJ = pseudochargeDensity_atom(V_PS,II,JJ,KK,xx(1),S);
		bJ_ref = pseudochargeDensity_atom(V_PS_ref,II,JJ,KK,xx(1),S);
		
		bJ = (-1/(4*pi))*bJ;
		bJ_ref = (-1/(4*pi))*bJ_ref;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
		% Calculate the gradient of pseudocharges
		% II = 1+2*S.FDn : size(V_PS,1)-2*S.FDn;
		% JJ = 1+2*S.FDn : size(V_PS,2)-2*S.FDn;
		% KK = 1+2*S.FDn : size(V_PS,3)-2*S.FDn;
		%[dbJ_x, dbJ_y, dbJ_z] = dpseudopot(bJ,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);
		%[dbJ_ref_x, dbJ_ref_y, dbJ_ref_z] = dpseudopot(bJ_ref,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);
		[dVJ_x, dVJ_y, dVJ_z] = dpseudopot(V_PS,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);
		[dVJ_ref_x, dVJ_ref_y, dVJ_ref_z] = dpseudopot(V_PS_ref,II,JJ,KK,XX_3D(II,JJ,KK),YY_3D(II,JJ,KK),ZZ_3D(II,JJ,KK),S);
		
		% Calculate local force and correction force components
		[II_rb,JJ_rb,KK_rb] = ndgrid(ii_s:ii_e,jj_s:jj_e,kk_s:kk_e);
		Rowcount_rb = (KK_rb-1)*S.Nx*S.Ny + (JJ_rb-1)*S.Nx + II_rb;

		% fl_x = sum(sum(sum(S.W(Rowcount_rb).*dbJ_x(II,JJ,KK).*(S.phi(Rowcount_rb)))));
		% fl_y = sum(sum(sum(S.W(Rowcount_rb).*dbJ_y(II,JJ,KK).*(S.phi(Rowcount_rb)))));
		% fl_z = sum(sum(sum(S.W(Rowcount_rb).*dbJ_z(II,JJ,KK).*(S.phi(Rowcount_rb)))));
		% 
		% fl_x = sum(sum(sum(S.W(Rowcount_rb).*dbJ_x(II,JJ,KK).*(S.phi(Rowcount_rb)-V_PS(II,JJ,KK)))));
		% fl_y = sum(sum(sum(S.W(Rowcount_rb).*dbJ_y(II,JJ,KK).*(S.phi(Rowcount_rb)-V_PS(II,JJ,KK)))));
		% fl_z = sum(sum(sum(S.W(Rowcount_rb).*dbJ_z(II,JJ,KK).*(S.phi(Rowcount_rb)-V_PS(II,JJ,KK)))));

		fl_x = - sum(sum(sum(S.W(Rowcount_rb).*bJ(II,JJ,KK).*(Dphi_x(Rowcount_rb) - dVJ_x(II,JJ,KK)))));
		fl_y = - sum(sum(sum(S.W(Rowcount_rb).*bJ(II,JJ,KK).*(Dphi_y(Rowcount_rb) - dVJ_y(II,JJ,KK)))));
		fl_z = - sum(sum(sum(S.W(Rowcount_rb).*bJ(II,JJ,KK).*(Dphi_z(Rowcount_rb) - dVJ_z(II,JJ,KK)))));
		
		% form 1:
		% fc_x = 0.5*sum(sum(sum(S.W(Rowcount_rb).*( dbJ_ref_x(II,JJ,KK) .* (S.V_c(Rowcount_rb)-V_PS_ref(II,JJ,KK)) +...%
		%     + dbJ_x(II,JJ,KK) .* (S.V_c(Rowcount_rb)+V_PS(II,JJ,KK)) +...%
		%     + (dVJ_ref_x(II,JJ,KK)-dVJ_x(II,JJ,KK)).*(S.b_ref(Rowcount_rb)+S.b(Rowcount_rb)) +...
		%     - bJ_ref(II,JJ,KK).*dVJ_ref_x(II,JJ,KK) + bJ(II,JJ,KK).*dVJ_x(II,JJ,KK)  ) )));%
		% fc_y = 0.5*sum(sum(sum(S.W(Rowcount_rb).*( dbJ_ref_y(II,JJ,KK) .* (S.V_c(Rowcount_rb)-V_PS_ref(II,JJ,KK)) +...%
		%     + dbJ_y(II,JJ,KK) .* (S.V_c(Rowcount_rb)+V_PS(II,JJ,KK)) +...%
		%     + (dVJ_ref_y(II,JJ,KK)-dVJ_y(II,JJ,KK)).*(S.b_ref(Rowcount_rb)+S.b(Rowcount_rb)) +...
		%     - bJ_ref(II,JJ,KK).*dVJ_ref_y(II,JJ,KK) + bJ(II,JJ,KK).*dVJ_y(II,JJ,KK)  ) )));%
		% fc_z = 0.5*sum(sum(sum(S.W(Rowcount_rb).*( dbJ_ref_z(II,JJ,KK) .* (S.V_c(Rowcount_rb)-V_PS_ref(II,JJ,KK)) +...%
		%     + dbJ_z(II,JJ,KK) .* (S.V_c(Rowcount_rb)+V_PS(II,JJ,KK)) +...%
		%     + (dVJ_ref_z(II,JJ,KK)-dVJ_z(II,JJ,KK)).*(S.b_ref(Rowcount_rb)+S.b(Rowcount_rb)) +...
		%      - bJ_ref(II,JJ,KK).*dVJ_ref_z(II,JJ,KK) + bJ(II,JJ,KK).*dVJ_z(II,JJ,KK) ) )));%
		% form 2:
		% fc_x = 0.5*sum(sum(sum(S.W(Rowcount_rb).*( (S.b(Rowcount_rb)+S.b_ref(Rowcount_rb)) .* (dVJ_ref_x(II,JJ,KK)-dVJ_x(II,JJ,KK)) +...
		%     + (dbJ_ref_x(II,JJ,KK) + dbJ_x(II,JJ,KK)) .* S.V_c(Rowcount_rb) ) )));
		% fc_y = 0.5*sum(sum(sum(S.W(Rowcount_rb).*( (S.b(Rowcount_rb)+S.b_ref(Rowcount_rb)) .* (dVJ_ref_y(II,JJ,KK)-dVJ_y(II,JJ,KK)) +...
		%     + (dbJ_ref_y(II,JJ,KK) + dbJ_y(II,JJ,KK)) .* S.V_c(Rowcount_rb) ) )));
		% fc_z = 0.5*sum(sum(sum(S.W(Rowcount_rb).*( (S.b(Rowcount_rb)+S.b_ref(Rowcount_rb)) .* (dVJ_ref_z(II,JJ,KK)-dVJ_z(II,JJ,KK)) +...
		%     + (dbJ_ref_z(II,JJ,KK) + dbJ_z(II,JJ,KK)) .* S.V_c(Rowcount_rb) ) )));
		% form 3:
		 % fc_x = 0.5*sum(sum(sum(S.W(Rowcount_rb).*( -(Dpseudo_x(Rowcount_rb)) .* (V_PS_ref(II,JJ,KK)-V_PS(II,JJ,KK)) +...
		 %    + (dbJ_ref_x(II,JJ,KK) + dbJ_x(II,JJ,KK)) .* S.V_c(Rowcount_rb) ) )));
		 % fc_y = 0.5*sum(sum(sum(S.W(Rowcount_rb).*( -(Dpseudo_y(Rowcount_rb)) .* (V_PS_ref(II,JJ,KK)-V_PS(II,JJ,KK)) +...
		 %    + (dbJ_ref_y(II,JJ,KK) + dbJ_y(II,JJ,KK)) .* S.V_c(Rowcount_rb) ) )));
		 % fc_z = 0.5*sum(sum(sum(S.W(Rowcount_rb).*( -(Dpseudo_z(Rowcount_rb)) .* (V_PS_ref(II,JJ,KK)-V_PS(II,JJ,KK)) +...
		 %    + (dbJ_ref_z(II,JJ,KK) + dbJ_z(II,JJ,KK)) .* S.V_c(Rowcount_rb) ) )));
		 % form 4:
		% fc_x = - 0.5*sum(sum(sum(S.W(Rowcount_rb).*( (Dpseudo_x(Rowcount_rb)) .* (V_PS_ref(II,JJ,KK)-V_PS(II,JJ,KK)) ...
		%      + (bJ(II,JJ,KK) + bJ_ref(II,JJ,KK)) .* DVc_x(Rowcount_rb) ) )));
		% fc_y = - 0.5*sum(sum(sum(S.W(Rowcount_rb).*( (Dpseudo_y(Rowcount_rb)) .* (V_PS_ref(II,JJ,KK)-V_PS(II,JJ,KK)) ...
		%      + (bJ(II,JJ,KK) + bJ_ref(II,JJ,KK)) .* DVc_y(Rowcount_rb) ) )));
		% fc_z = - 0.5*sum(sum(sum(S.W(Rowcount_rb).*( (Dpseudo_z(Rowcount_rb)) .* (V_PS_ref(II,JJ,KK)-V_PS(II,JJ,KK)) ...
		%      + (bJ(II,JJ,KK) + bJ_ref(II,JJ,KK)) .* DVc_z(Rowcount_rb) ) )));
		 % form 5:
		fc_x = 0.5*sum(sum(sum(S.W(Rowcount_rb).*( (S.b(Rowcount_rb)+S.b_ref(Rowcount_rb)) .* (dVJ_ref_x(II,JJ,KK)-dVJ_x(II,JJ,KK)) +...
			 - (bJ(II,JJ,KK) + bJ_ref(II,JJ,KK)) .* DVc_x(Rowcount_rb) ) )));
		fc_y = 0.5*sum(sum(sum(S.W(Rowcount_rb).*( (S.b(Rowcount_rb)+S.b_ref(Rowcount_rb)) .* (dVJ_ref_y(II,JJ,KK)-dVJ_y(II,JJ,KK)) +...
			 - (bJ(II,JJ,KK) + bJ_ref(II,JJ,KK)) .* DVc_y(Rowcount_rb) ) )));
		fc_z = 0.5*sum(sum(sum(S.W(Rowcount_rb).*( (S.b(Rowcount_rb)+S.b_ref(Rowcount_rb)) .* (dVJ_ref_z(II,JJ,KK)-dVJ_z(II,JJ,KK)) +...
			 - (bJ(II,JJ,KK) + bJ_ref(II,JJ,KK)) .* DVc_z(Rowcount_rb) ) )));
		

		% Perform Rotations for cychel
		if(S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
			fac1 = (y0-y0_i)/S.L2;
			fac2 = (z0-z0_i)/S.L3;
			ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
			P = ROT*[fl_x;fl_y;fl_z];
			fl_x = P(1); fl_y = P(2); fl_z = P(3);
			P = ROT*[fc_x;fc_y;fc_z];
			fc_x = P(1); fc_y = P(2); fc_z = P(3);
		end


		force_local(JJ_a,1) = force_local(JJ_a,1) + fl_x;
		force_local(JJ_a,2) = force_local(JJ_a,2) + fl_y;
		force_local(JJ_a,3) = force_local(JJ_a,3) + fl_z;

		force_corr(JJ_a,1) = force_corr(JJ_a,1) + fc_x;
		force_corr(JJ_a,2) = force_corr(JJ_a,2) + fc_y;
		force_corr(JJ_a,3) = force_corr(JJ_a,3) + fc_z;
	end
	
	% Check if same type of atoms are over
	if count_typ_atms == S.Atm(count_typ).n_atm_typ
		count_typ_atms = 1;
		count_typ = count_typ + 1;
	else
		count_typ_atms = count_typ_atms + 1;
	end

end % end of loop over atoms


fprintf(' local force calculation: %.3f s\n', toc(tic_locforces));

%**********************************************************************
%*                   Calculate nonlocal atomic force                  *
%**********************************************************************

for kpt = 1:S.tnkpt
    if (kpt(1) == 0 && kpt(2) == 0 && kpt(3) == 0)
        fac = 1.0;
    else
        fac = 1.0i;
    end
    kpt_vec = S.kptgrid(kpt,:);

    % Calculate gradient of psi   
    Dpsi_x = zeros(S.nspinor *S.N,S.Nev);
    Dpsi_y = zeros(S.nspinor *S.N,S.Nev);
    Dpsi_z = zeros(S.nspinor *S.N,S.Nev);
    for spinor = 1:S.nspinor
        ndrange = (1+(spinor-1)*S.N:spinor*S.N); 
        Dpsi_x(ndrange,:) = blochGradient(S,kpt_vec,1)*S.psi(ndrange,:,kpt);
        Dpsi_y(ndrange,:) = blochGradient(S,kpt_vec,2)*S.psi(ndrange,:,kpt);
        Dpsi_z(ndrange,:) = blochGradient(S,kpt_vec,3)*S.psi(ndrange,:,kpt);
    end

    for spinor = 1:S.nspinor
        sigma = (-1)^(spinor-1);
        shift = (spinor-1)*S.N;         % for selecting each spinor, spinor=1 shift = 0, spinor=2, shift = S.N
        shift2 = (2-spinor)*S.N;        % for selecting the other spin channel, spinor=1 shift2 = S.N, spinor=2,shift2=0  
        nsshift = (spinor-1)*S.tnkpt*(S.spin_typ == 1);

        % scalar relativistic part
        force_atm = zeros(S.n_atm,3);
        for JJ_a = 1:S.n_atm % loop over all atoms
            % Calculate nonlocal components of the force acting on atom JJ_a
            integral_1 = zeros(S.Atom(JJ_a).angnum,S.Nev);
            integral_2_x = zeros(S.Atom(JJ_a).angnum,S.Nev);
            integral_2_y = zeros(S.Atom(JJ_a).angnum,S.Nev);
            integral_2_z = zeros(S.Atom(JJ_a).angnum,S.Nev);
            for img = 1:S.Atom(JJ_a).n_image_rc
                img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
                phase_fac = exp(dot(kpt_vec,img_disp*fac));
                ChiW = transpose(bsxfun(@times, conj(S.Atom(JJ_a).rcImage(img).Chi_mat), S.W(S.Atom(JJ_a).rcImage(img).rc_pos)));
                integral_1 = integral_1 + conj(ChiW) * conj(S.psi(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:,kpt)) * conj(phase_fac);
                integral_2_x = integral_2_x + ChiW * (Dpsi_x(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)) * (phase_fac);
                integral_2_y = integral_2_y + ChiW * (Dpsi_y(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)) * (phase_fac);
                integral_2_z = integral_2_z + ChiW * (Dpsi_z(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)) * (phase_fac);
            end
            tf_x = transpose(S.Atom(JJ_a).gamma_Jl) * real(integral_1.*integral_2_x) * S.occ(:,kpt+nsshift);
            tf_y = transpose(S.Atom(JJ_a).gamma_Jl) * real(integral_1.*integral_2_y) * S.occ(:,kpt+nsshift);
            tf_z = transpose(S.Atom(JJ_a).gamma_Jl) * real(integral_1.*integral_2_z) * S.occ(:,kpt+nsshift);
            force_atm(JJ_a,:) = force_atm(JJ_a,:) + [tf_x tf_y tf_z];
        end
        force_nloc = force_nloc - S.occfac*2*S.wkpt(kpt)*force_atm;
        
        % spin-orbit coupling part 1
        force_atm = zeros(S.n_atm,3);
        for JJ_a = 1:S.n_atm % loop over all atoms
            if S.Atm(S.Atom(JJ_a).count_typ).pspsoc == 0
                continue;
            end
            ncol_term1 = S.Atom(JJ_a).ncol_term1;
            soindx = S.Atom(JJ_a).term1_index_so(1:ncol_term1);
            
            % Calculate nonlocal components of the force acting on atom JJ_a
            integral_1 = zeros(ncol_term1,S.Nev);
            integral_2_x = zeros(ncol_term1,S.Nev);
            integral_2_y = zeros(ncol_term1,S.Nev);
            integral_2_z = zeros(ncol_term1,S.Nev);
            for img = 1:S.Atom(JJ_a).n_image_rc
                img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
                phase_fac = exp(dot(kpt_vec,img_disp*fac));
                ChiW = transpose(bsxfun(@times, conj(S.Atom(JJ_a).rcImage(img).Chiso_mat(:,soindx)), S.W(S.Atom(JJ_a).rcImage(img).rc_pos)));
                integral_1 = integral_1 + conj(ChiW) * conj(S.psi(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:,kpt)) * conj(phase_fac);
                integral_2_x = integral_2_x + ChiW * (Dpsi_x(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)) * (phase_fac);
                integral_2_y = integral_2_y + ChiW * (Dpsi_y(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)) * (phase_fac);
                integral_2_z = integral_2_z + ChiW * (Dpsi_z(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)) * (phase_fac);
            end
            tf_x = transpose(sigma*S.Atom(JJ_a).term1_gammaso_Jl(1:ncol_term1)) * real(integral_1.*integral_2_x) * S.occ(:,kpt);
            tf_y = transpose(sigma*S.Atom(JJ_a).term1_gammaso_Jl(1:ncol_term1)) * real(integral_1.*integral_2_y) * S.occ(:,kpt);
            tf_z = transpose(sigma*S.Atom(JJ_a).term1_gammaso_Jl(1:ncol_term1)) * real(integral_1.*integral_2_z) * S.occ(:,kpt);
            force_atm(JJ_a,:) = force_atm(JJ_a,:) + [tf_x tf_y tf_z];
        end
        force_nloc = force_nloc - S.occfac*2*S.wkpt(kpt)*force_atm;
        
        % spin-orbit coupling part 2
        force_atm = zeros(S.n_atm,3);
        for JJ_a = 1:S.n_atm % loop over all atoms
            if S.Atm(S.Atom(JJ_a).count_typ).pspsoc == 0
                continue;
            end
            ncol_term2 = S.Atom(JJ_a).ncol_term2;
            if spinor == 1
                soindx1 = S.Atom(JJ_a).term2_index_so(1:ncol_term2)+1;
                soindx2 = S.Atom(JJ_a).term2_index_so(1:ncol_term2);
            else 
                soindx1 = S.Atom(JJ_a).term2_index_so(1:ncol_term2);
                soindx2 = S.Atom(JJ_a).term2_index_so(1:ncol_term2)+1;
            end
            
            % Calculate nonlocal components of the force acting on atom JJ_a
            integral_1 = zeros(ncol_term2,S.Nev);
            integral_2_x = zeros(ncol_term2,S.Nev);
            integral_2_y = zeros(ncol_term2,S.Nev);
            integral_2_z = zeros(ncol_term2,S.Nev);
            
            for img = 1:S.Atom(JJ_a).n_image_rc
                img_disp = S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates;
                phase_fac = exp(dot(kpt_vec,img_disp*fac));
                ChiW1 = transpose(bsxfun(@times, conj(S.Atom(JJ_a).rcImage(img).Chiso_mat(:,soindx1)), S.W(S.Atom(JJ_a).rcImage(img).rc_pos)));
                ChiW2 = transpose(bsxfun(@times, S.Atom(JJ_a).rcImage(img).Chiso_mat(:,soindx2), S.W(S.Atom(JJ_a).rcImage(img).rc_pos)));
                
                integral_1 = integral_1 + ChiW2 * conj(S.psi(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:,kpt)) * conj(phase_fac);
                integral_2_x = integral_2_x + ChiW1 * (Dpsi_x(S.Atom(JJ_a).rcImage(img).rc_pos+shift2,:)) * (phase_fac);
                integral_2_y = integral_2_y + ChiW1 * (Dpsi_y(S.Atom(JJ_a).rcImage(img).rc_pos+shift2,:)) * (phase_fac);
                integral_2_z = integral_2_z + ChiW1 * (Dpsi_z(S.Atom(JJ_a).rcImage(img).rc_pos+shift2,:)) * (phase_fac);
            end
            tf_x = transpose(S.Atom(JJ_a).term2_gammaso_Jl(1:ncol_term2)) * real(integral_1.*integral_2_x) * S.occ(:,kpt);
            tf_y = transpose(S.Atom(JJ_a).term2_gammaso_Jl(1:ncol_term2)) * real(integral_1.*integral_2_y) * S.occ(:,kpt);
            tf_z = transpose(S.Atom(JJ_a).term2_gammaso_Jl(1:ncol_term2)) * real(integral_1.*integral_2_z) * S.occ(:,kpt);
            force_atm(JJ_a,:) = force_atm(JJ_a,:) + [tf_x tf_y tf_z];
        end
        force_nloc = force_nloc - S.occfac*2*S.wkpt(kpt)*force_atm;
    end
end


% Total force

force = force_local + force_corr + force_nloc + force_xc;
% force_local
% force_corr
% force_nloc
if S.cell_typ == 2
	force = force*S.grad_T; % Convert forces from lattice to cartesian coordinates
end

force = real(force);

if S.d3Flag == 1 
	force = force - S.d3gradAtomPos;
end
end



