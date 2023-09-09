function pressure = evaluatePressure(S)
% @brief    Function to calculate pressure in periodic systems (O(N^3))
% @authors
%          Abhiraj Sharma <asharma424@gatech.edu>
%          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @references
%             "On the calculation of the stress tensor in real-space Kohn-Sham
%              density functional theory (Sharma et. al. 2018)"
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%=================================================================================

% Component from Energy related terms calculated explicitly
Dphi_x = S.grad_1*(S.phi);
Dphi_y = S.grad_2*(S.phi);
Dphi_z = S.grad_3*(S.phi);

P_eng = 0.0;
ks = 1;
for spin = 1:S.nspin
	for kpt = 1:S.tnkpt
		P_eng = P_eng - S.occfac*2*S.wkpt(kpt)*sum(S.EigVal(:,ks).*S.occ(:,ks)) ;
		ks = ks + 1;
	end
end

if S.spin_typ == 0
	Drho_x = S.grad_1 * (S.rho+S.rho_Tilde_at);
	Drho_y = S.grad_2 * (S.rho+S.rho_Tilde_at);
	Drho_z = S.grad_3 * (S.rho+S.rho_Tilde_at);
else
	rho = S.rho;
	rho(:,1) = rho(:,1)+S.rho_Tilde_at;
	rho(:,2) = rho(:,2)+0.5*S.rho_Tilde_at;
	rho(:,3) = rho(:,3)+0.5*S.rho_Tilde_at;
	Drho_x = S.grad_1*rho;
	Drho_y = S.grad_2*rho;
	Drho_z = S.grad_3*rho;
end

P_eng = P_eng + 0.5 * sum(S.rho(:,1) .* S.phi .* S.W) + ...
	(1.5 * sum(S.b .* S.phi .* S.W))  + (3*S.E_corr) - 3*S.Eself + ...
	(1/4/pi) * sum( S.W .* (S.lapc_T(1,1)*Dphi_x.*Dphi_x + S.lapc_T(2,2)*Dphi_y.*Dphi_y + S.lapc_T(3,3)*Dphi_z.*Dphi_z +...
							S.lapc_T(1,2)*Dphi_x.*Dphi_y + S.lapc_T(2,3)*Dphi_y.*Dphi_z + S.lapc_T(1,3)*Dphi_x.*Dphi_z )) ;

P_nlcc= 0;
if S.NLCC_flag
    S_T = transpose(S.lat_uvec);
    count_typ = 1;
    count_typ_atms = 1;
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

            drho_Tilde_at_1 = S.grad_T(1,1)*drho_Tilde_at_x + S.grad_T(2,1)*drho_Tilde_at_y + S.grad_T(3,1)*drho_Tilde_at_z;
            drho_Tilde_at_2 = S.grad_T(1,2)*drho_Tilde_at_x + S.grad_T(2,2)*drho_Tilde_at_y + S.grad_T(3,2)*drho_Tilde_at_z;
            drho_Tilde_at_3 = S.grad_T(1,3)*drho_Tilde_at_x + S.grad_T(2,3)*drho_Tilde_at_y + S.grad_T(3,3)*drho_Tilde_at_z;

            % Calculate xc 2nd term stress components
            [II_rb,JJ_rb,KK_rb] = ndgrid(ii_s:ii_e,jj_s:jj_e,kk_s:kk_e);
            Rowcount_rb = (KK_rb-1)*S.Nx*S.Ny + (JJ_rb-1)*S.Nx + II_rb;
            [xr,yr,zr] = ndgrid((ii_s-1:ii_e-1)*S.dx - x0_i,(jj_s-1:jj_e-1)*S.dy - y0_i,(kk_s-1:kk_e-1)*S.dz - z0_i) ;
            x1 = S_T(1,1)*xr + S_T(1,2)*yr + S_T(1,3)*zr;
            y1 = S_T(2,1)*xr + S_T(2,2)*yr + S_T(2,3)*zr;
            z1 = S_T(3,1)*xr + S_T(3,2)*yr + S_T(3,3)*zr;
            if S.spin_typ == 0
                P_nlcc = P_nlcc + sum(sum(sum( drho_Tilde_at_1(II,JJ,KK) .* x1 .* ( S.Vxc(Rowcount_rb) ) .* S.W(Rowcount_rb) )));
                P_nlcc = P_nlcc + sum(sum(sum( drho_Tilde_at_2(II,JJ,KK) .* y1 .* ( S.Vxc(Rowcount_rb) ) .* S.W(Rowcount_rb) )));
                P_nlcc = P_nlcc + sum(sum(sum( drho_Tilde_at_3(II,JJ,KK) .* z1 .* ( S.Vxc(Rowcount_rb) ) .* S.W(Rowcount_rb) )));
            else
                vxc = S.Vxc(:,1)+S.Vxc(:,2);
                P_nlcc = P_nlcc + sum(sum(sum( 0.5*drho_Tilde_at_1(II,JJ,KK) .* x1 .* ( vxc(Rowcount_rb)  ) .* S.W(Rowcount_rb) )));
                P_nlcc = P_nlcc + sum(sum(sum( 0.5*drho_Tilde_at_2(II,JJ,KK) .* y1 .* (  vxc(Rowcount_rb)  ) .* S.W(Rowcount_rb) )));
                P_nlcc = P_nlcc + sum(sum(sum( 0.5*drho_Tilde_at_3(II,JJ,KK) .* z1 .* (  vxc(Rowcount_rb)  ) .* S.W(Rowcount_rb) )));
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

% Contribution from exchange-correlation
if S.spin_typ == 0
	P_eng = P_eng - sum(S.rho .* S.Vxc .* S.W) + (3*S.Exc) + ...
			- S.W' * (S.dvxcdgrho .* (S.lapc_T(1,1)*Drho_x.*Drho_x + S.lapc_T(2,2)*Drho_y.*Drho_y + S.lapc_T(3,3)*Drho_z.*Drho_z +...
										 S.lapc_T(1,2)*Drho_x.*Drho_y + S.lapc_T(2,3)*Drho_y.*Drho_z + S.lapc_T(1,3)*Drho_z.*Drho_x ));
    if (S.countPotential > 0) &&(S.xc == 4) % add metaGGA pressure term
		P_eng = P_eng - sum(S.W .* S.VxcScan3 .* S.tau);
        for kpt = 1:S.tnkpt
			kpt_vec = S.kptgrid(kpt,:);

		    Dpsi_x = blochGradient(S,kpt_vec,1)*S.psi(:,:,kpt);
		    Dpsi_y = blochGradient(S,kpt_vec,2)*S.psi(:,:,kpt);
		    Dpsi_z = blochGradient(S,kpt_vec,3)*S.psi(:,:,kpt);

		    TDpsi_1 = S.grad_T(1,1)*Dpsi_x + S.grad_T(2,1)*Dpsi_y + S.grad_T(3,1)*Dpsi_z;
		    TDpsi_2 = S.grad_T(1,2)*Dpsi_x + S.grad_T(2,2)*Dpsi_y + S.grad_T(3,2)*Dpsi_z;
		    TDpsi_3 = S.grad_T(1,3)*Dpsi_x + S.grad_T(2,3)*Dpsi_y + S.grad_T(3,3)*Dpsi_z;
		    TDcpsi_1 = conj(TDpsi_1);
		    TDcpsi_2 = conj(TDpsi_2);
		    TDcpsi_3 = conj(TDpsi_3);

			P_eng = P_eng - real(S.occfac * S.wkpt(kpt) * (transpose(S.W) * ... 
				(S.VxcScan3.*TDcpsi_1.*TDpsi_1)) * S.occ(:,kpt));
			P_eng = P_eng - real(S.occfac * S.wkpt(kpt) * (transpose(S.W) * ... 
				(S.VxcScan3.*TDcpsi_2.*TDpsi_2)) * S.occ(:,kpt));
			P_eng = P_eng - real(S.occfac * S.wkpt(kpt) * (transpose(S.W) * ... 
				(S.VxcScan3.*TDcpsi_3.*TDpsi_3)) * S.occ(:,kpt));
        end
    end
else
    rho = S.rho;
	P_eng = P_eng - sum(sum(S.Vxc.*rho(:,2:3),2).*S.W) + (3*S.Exc) + ...
			- sum(S.W' * (S.dvxcdgrho .* (S.lapc_T(1,1)*Drho_x.*Drho_x + S.lapc_T(2,2)*Drho_y.*Drho_y + S.lapc_T(3,3)*Drho_z.*Drho_z +...
											 S.lapc_T(1,2)*Drho_x.*Drho_y + S.lapc_T(2,3)*Drho_y.*Drho_z + S.lapc_T(1,3)*Drho_z.*Drho_x )));
end
% Components from gradient terms

% initialization

P_elec = 0;
P_corr = 0;
%psdfilepath = sprintf('%s/PseudopotFiles',S.inputfile_path);
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
		
		% Indices of closest grid point to atom
		pos_ii = round((x0_i-S.xin) / S.dx) + 1;
		pos_jj = round((y0_i-S.yin) / S.dy) + 1;
		pos_kk = round((z0_i-S.zin) / S.dz) + 1;
		
		% Starting and ending indices of b-region
		ii_s = pos_ii - ceil(S.Atm(count_typ).rb_x/S.dx+0.5);
		ii_e = pos_ii + ceil(S.Atm(count_typ).rb_x/S.dx+0.5);
		jj_s = pos_jj - ceil(S.Atm(count_typ).rb_y/S.dy+0.5);
		jj_e = pos_jj + ceil(S.Atm(count_typ).rb_y/S.dy+0.5);
		kk_s = pos_kk - ceil(S.Atm(count_typ).rb_z/S.dz+0.5);
		kk_e = pos_kk + ceil(S.Atm(count_typ).rb_z/S.dz+0.5);
		
		% Check if the b-region is inside the domain in Dirichlet BC
		% direction
		isInside = (S.BCx == 0 || (S.BCx == 1 && (ii_s>1) && (ii_e<S.Nx))) && ...
		   (S.BCy == 0 || (S.BCy == 1 && (jj_s>1) && (jj_e<S.Ny))) && ...
		   (S.BCz == 0 || (S.BCz == 1 && (kk_s>1) && (kk_e<S.Nz)));
		assert(isInside,'ERROR: Atom too close to boundary for b calculation');
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
		V_PS = zeros(size(dd));
		IsLargeThanRmax = dd > S.Atm(count_typ).r_grid_vloc(end);
		V_PS(IsLargeThanRmax) = -S.Atm(count_typ).Z;
		V_PS(~IsLargeThanRmax) = interp1(S.Atm(count_typ).r_grid_vloc, S.Atm(count_typ).r_grid_vloc.*S.Atm(count_typ).Vloc, dd(~IsLargeThanRmax), 'spline');

		V_PS = V_PS./dd;
		V_PS(dd<S.Atm(count_typ).r_grid_vloc(2)) = S.Atm(count_typ).Vloc(1); % WARNING
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
		
		% Calculate the gradient of pseudocharges
		dbJ_x = zeros(size(V_PS)); dbJ_y = zeros(size(V_PS)); dbJ_z = zeros(size(V_PS));
		dbJ_ref_x = zeros(size(V_PS)); dbJ_ref_y = zeros(size(V_PS)); dbJ_ref_z = zeros(size(V_PS));
		dVJ_x = zeros(size(V_PS)); dVJ_y = zeros(size(V_PS)); dVJ_z = zeros(size(V_PS));
		dVJ_ref_x = zeros(size(V_PS)); dVJ_ref_y = zeros(size(V_PS)); dVJ_ref_z = zeros(size(V_PS));
				
		II = 1+2*S.FDn : size(V_PS,1)-2*S.FDn;
		JJ = 1+2*S.FDn : size(V_PS,2)-2*S.FDn;
		KK = 1+2*S.FDn : size(V_PS,3)-2*S.FDn;
		
		for p = 1:S.FDn
			dbJ_x(II,JJ,KK) = dbJ_x(II,JJ,KK) + S.w1(p+1)/S.dx*(bJ(II+p,JJ,KK)-bJ(II-p,JJ,KK));
			dbJ_y(II,JJ,KK) = dbJ_y(II,JJ,KK) + S.w1(p+1)/S.dy*(bJ(II,JJ+p,KK)-bJ(II,JJ-p,KK));
			dbJ_z(II,JJ,KK) = dbJ_z(II,JJ,KK) + S.w1(p+1)/S.dz*(bJ(II,JJ,KK+p)-bJ(II,JJ,KK-p));
			dbJ_ref_x(II,JJ,KK) = dbJ_ref_x(II,JJ,KK) + S.w1(p+1)/S.dx*(bJ_ref(II+p,JJ,KK)-bJ_ref(II-p,JJ,KK));
			dbJ_ref_y(II,JJ,KK) = dbJ_ref_y(II,JJ,KK) + S.w1(p+1)/S.dy*(bJ_ref(II,JJ+p,KK)-bJ_ref(II,JJ-p,KK));
			dbJ_ref_z(II,JJ,KK) = dbJ_ref_z(II,JJ,KK) + S.w1(p+1)/S.dz*(bJ_ref(II,JJ,KK+p)-bJ_ref(II,JJ,KK-p));
			dVJ_x(II,JJ,KK) = dVJ_x(II,JJ,KK) + S.w1(p+1)/S.dx*(V_PS(II+p,JJ,KK)-V_PS(II-p,JJ,KK));
			dVJ_y(II,JJ,KK) = dVJ_y(II,JJ,KK) + S.w1(p+1)/S.dy*(V_PS(II,JJ+p,KK)-V_PS(II,JJ-p,KK));
			dVJ_z(II,JJ,KK) = dVJ_z(II,JJ,KK) + S.w1(p+1)/S.dz*(V_PS(II,JJ,KK+p)-V_PS(II,JJ,KK-p));
			dVJ_ref_x(II,JJ,KK) = dVJ_ref_x(II,JJ,KK) + S.w1(p+1)/S.dx*(V_PS_ref(II+p,JJ,KK)-V_PS_ref(II-p,JJ,KK));
			dVJ_ref_y(II,JJ,KK) = dVJ_ref_y(II,JJ,KK) + S.w1(p+1)/S.dy*(V_PS_ref(II,JJ+p,KK)-V_PS_ref(II,JJ-p,KK));
			dVJ_ref_z(II,JJ,KK) = dVJ_ref_z(II,JJ,KK) + S.w1(p+1)/S.dz*(V_PS_ref(II,JJ,KK+p)-V_PS_ref(II,JJ,KK-p));
		end
		
		% Calculate local (electrostatics) pressure components:
		[II_rb,JJ_rb,KK_rb] = ndgrid(ii_s:ii_e,jj_s:jj_e,kk_s:kk_e);
		Rowcount_rb = (KK_rb-1)*S.Nx*S.Ny + (JJ_rb-1)*S.Nx + II_rb;
		[xr,yr,zr] = ndgrid((ii_s-1:ii_e-1)*S.dx - x0_i,(jj_s-1:jj_e-1)*S.dy - y0_i,(kk_s-1:kk_e-1)*S.dz - z0_i) ;
		
		P_elec = P_elec + sum(sum(sum( dVJ_x(II,JJ,KK) .* xr .* ( - 0.5 * bJ(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
		P_elec = P_elec + sum(sum(sum( dVJ_y(II,JJ,KK) .* yr .* ( - 0.5 * bJ(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
		P_elec = P_elec + sum(sum(sum( dVJ_z(II,JJ,KK) .* zr .* ( - 0.5 * bJ(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
		
		%         P_Eself = P_Eself + sum(sum(sum( dVJ_x(II,JJ,KK) .* xr .* ( 0.5 * bJ(II,JJ,KK) ) .* S.W(Rowcount_rb) ))) + sum(sum(sum( dbJ_x(II,JJ,KK) .* xr .* (0.5 * V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
		%         P_Eself = P_Eself + sum(sum(sum( dVJ_y(II,JJ,KK) .* yr .* ( 0.5 * bJ(II,JJ,KK) ) .* S.W(Rowcount_rb) ))) + sum(sum(sum( dbJ_y(II,JJ,KK) .* yr .* (0.5 * V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
		%         P_Eself = P_Eself + sum(sum(sum( dVJ_z(II,JJ,KK) .* zr .* ( 0.5 * bJ(II,JJ,KK) ) .* S.W(Rowcount_rb) ))) + sum(sum(sum( dbJ_z(II,JJ,KK) .* zr .* (0.5 * V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
		%
		%         P_Eself = P_Eself + 0.5 * sum(sum(sum( Dself_x(II,JJ,KK) .*xr .* S.W(Rowcount_rb))));
		%         P_Eself = P_Eself + 0.5 * sum(sum(sum( Dself_y(II,JJ,KK) .*yr .* S.W(Rowcount_rb))));
		%         P_Eself = P_Eself + 0.5 * sum(sum(sum( Dself_z(II,JJ,KK) .*zr .* S.W(Rowcount_rb))));
		
		P_elec = P_elec + sum(sum(sum( dbJ_x(II,JJ,KK) .* xr .* ( S.phi(Rowcount_rb) -0.5 * V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
		P_elec = P_elec + sum(sum(sum( dbJ_y(II,JJ,KK) .* yr .* ( S.phi(Rowcount_rb) -0.5 * V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
		P_elec = P_elec + sum(sum(sum( dbJ_z(II,JJ,KK) .* zr .* ( S.phi(Rowcount_rb) -0.5 * V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
		
		%P_elec = P_elec + sum(sum(sum( dbJ_x(II,JJ,KK) .* xr .* ( S.phi(Rowcount_rb)) .* S.W(Rowcount_rb) )));
		%P_elec = P_elec + sum(sum(sum( dbJ_y(II,JJ,KK) .* yr .* ( S.phi(Rowcount_rb)) .* S.W(Rowcount_rb) )));
		%P_elec = P_elec + sum(sum(sum( dbJ_z(II,JJ,KK) .* zr .* ( S.phi(Rowcount_rb)) .* S.W(Rowcount_rb) )));
		
		%         P_elec = P_elec + sum(sum(sum( (Dbphi_x(II,JJ,KK) .* xp - dbJ_x(II,JJ,KK).*S.phi(Rowcount_rb)*x0_i).* S.W(Rowcount_rb)  )));
		%         P_elec = P_elec + sum(sum(sum( (Dbphi_y(II,JJ,KK) .* yp - dbJ_y(II,JJ,KK).*S.phi(Rowcount_rb)*y0_i).* S.W(Rowcount_rb) )));
		%         P_elec = P_elec + sum(sum(sum( (Dbphi_z(II,JJ,KK) .* zp - dbJ_z(II,JJ,KK).*S.phi(Rowcount_rb)*z0_i).* S.W(Rowcount_rb) )));
		
		P_corr = P_corr + 0.5 * sum(sum(sum( ( dbJ_x(II,JJ,KK) .* ( S.V_c(Rowcount_rb) + V_PS(II,JJ,KK) ) + ...
			dbJ_ref_x(II,JJ,KK) .* ( S.V_c(Rowcount_rb) - V_PS_ref(II,JJ,KK) ) + ...
			dVJ_ref_x(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ_ref(II,JJ,KK)) - ...
			dVJ_x(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ(II,JJ,KK)) ) .* xr .* S.W(Rowcount_rb) )));
		
		P_corr = P_corr + 0.5 * sum(sum(sum( ( dbJ_y(II,JJ,KK) .* ( S.V_c(Rowcount_rb) + V_PS(II,JJ,KK) ) + ...
			dbJ_ref_y(II,JJ,KK) .* ( S.V_c(Rowcount_rb) - V_PS_ref(II,JJ,KK) ) + ...
			dVJ_ref_y(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ_ref(II,JJ,KK)) - ...
			dVJ_y(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ(II,JJ,KK)) ) .* yr .* S.W(Rowcount_rb) )));
		
		P_corr = P_corr + 0.5 * sum(sum(sum( ( dbJ_z(II,JJ,KK) .* ( S.V_c(Rowcount_rb) + V_PS(II,JJ,KK) ) + ...
			dbJ_ref_z(II,JJ,KK) .* ( S.V_c(Rowcount_rb) - V_PS_ref(II,JJ,KK) ) + ...
			dVJ_ref_z(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ_ref(II,JJ,KK)) - ...
			dVJ_z(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ(II,JJ,KK)) ) .* zr .* S.W(Rowcount_rb) )));
		
		
		%         P_elec = P_elec + sum(sum(sum( dbJ_x(II,JJ,KK) .* xr .* ( S.phi(Rowcount_rb) - V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
		%         P_elec = P_elec + sum(sum(sum( dbJ_y(II,JJ,KK) .* yr .* ( S.phi(Rowcount_rb) - V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
		%         P_elec = P_elec + sum(sum(sum( dbJ_z(II,JJ,KK) .* zr .* ( S.phi(Rowcount_rb) - V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
		%
		%        P_corr = P_corr + 0.5 * sum(sum(sum( ( (dbJ_x(II,JJ,KK) + dbJ_ref_x(II,JJ,KK)).* S.V_c(Rowcount_rb) + ...
		%            (dVJ_ref_x(II,JJ,KK) - dVJ_x(II,JJ,KK)) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb)) ) .* xr .* S.W(Rowcount_rb) )));
		
		%         P_corr = P_corr + 0.5 * sum(sum(sum( ( (dbJ_y(II,JJ,KK) + dbJ_ref_y(II,JJ,KK)).* S.V_c(Rowcount_rb) + ...
		%            (dVJ_ref_y(II,JJ,KK) - dVJ_y(II,JJ,KK)) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb)) ) .* yr .* S.W(Rowcount_rb) )));
		
		%         P_corr = P_corr + 0.5 * sum(sum(sum( ( (dbJ_z(II,JJ,KK) + dbJ_ref_z(II,JJ,KK)).* S.V_c(Rowcount_rb) + ...
		%             (dVJ_ref_z(II,JJ,KK) - dVJ_z(II,JJ,KK)) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb)) ) .* zr .* S.W(Rowcount_rb) )));
		
	end
	
	
	% Check if same type of atoms are over
	if count_typ_atms == S.Atm(count_typ).n_atm_typ
		count_typ_atms = 1;
		count_typ = count_typ + 1;
	else
		count_typ_atms = count_typ_atms + 1;
	end
	
end % end of loop over atoms


%**********************************************************************
%*                   Calculate nonlocal Pressure                  *
%**********************************************************************
%Type-I (gradient on psi rather than on Chi)

P_nl = 0;    
for kpt = 1:S.tnkpt
    fac = 1.0i;
    kpt_vec = S.kptgrid(kpt,:);
    Dpsi_x = zeros(S.nspinor *S.N,S.Nev);
    Dpsi_y = zeros(S.nspinor *S.N,S.Nev);
    Dpsi_z = zeros(S.nspinor *S.N,S.Nev);
    for spinor = 1:S.nspinor
        ndrange = (1+(spinor-1)*S.N:spinor*S.N); 
        Dpsi_x(ndrange,:) = blochGradient(S,kpt_vec,1)*S.psi(ndrange,:,kpt);
        Dpsi_y(ndrange,:) = blochGradient(S,kpt_vec,2)*S.psi(ndrange,:,kpt);
        Dpsi_z(ndrange,:) = blochGradient(S,kpt_vec,3)*S.psi(ndrange,:,kpt);
    end

    for JJ_a = 1:S.n_atm % loop over all atoms
        for spinor = 1:S.nspinor
            sigma = (-1)^(spinor-1);
            shift = (spinor-1)*S.N;         % for selecting each spinor, spinor=1 shift = 0, spinor=2, shift = S.N
            nsshift = (spinor-1)*S.tnkpt*(S.spin_typ == 1);
            
            % Nonlocal scalar relativisic energy
            Chi_X_mult1 = zeros(S.Atom(JJ_a).angnum,S.Nev);
            for img = 1:S.Atom(JJ_a).n_image_rc
                phase_fac = (exp(dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)*fac)));
                Chi_X_mult1 = Chi_X_mult1 + transpose(bsxfun(@times, conj(S.Atom(JJ_a).rcImage(img).Chi_mat), S.W(S.Atom(JJ_a).rcImage(img).rc_pos))) * S.psi(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:,kpt) * phase_fac ;
            end
            
            P_nl = P_nl - S.occfac * S.wkpt(kpt) * transpose(S.Atom(JJ_a).gamma_Jl) * (Chi_X_mult1.*conj(Chi_X_mult1)) * S.occ(:,kpt+nsshift);
            
            % Nonlocal spin-orbit coupling term1 energy
            if S.Atm(S.Atom(JJ_a).count_typ).pspsoc == 1
                ncol_term1 = S.Atom(JJ_a).ncol_term1;
                soindx = S.Atom(JJ_a).term1_index_so(1:ncol_term1);
                Chiso_X_mult1 = zeros(ncol_term1,S.Nev);
                for img = 1:S.Atom(JJ_a).n_image_rc
                    phase_fac = (exp(dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)*fac)));
                    Chiso_X_mult1 = Chiso_X_mult1 + transpose(bsxfun(@times, conj(S.Atom(JJ_a).rcImage(img).Chiso_mat(:,soindx)), S.W(S.Atom(JJ_a).rcImage(img).rc_pos))) * S.psi(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:,kpt) * phase_fac ;
                end

                P_nl = P_nl - S.occfac * S.wkpt(kpt) * transpose(sigma*S.Atom(JJ_a).term1_gammaso_Jl(1:ncol_term1)) * (Chiso_X_mult1.*conj(Chiso_X_mult1)) * S.occ(:,kpt) ;
            end
        end
        
        % Nonlocal spin-orbit coupling term2 energy
        if S.Atm(S.Atom(JJ_a).count_typ).pspsoc == 1
            ncol_term2 = S.Atom(JJ_a).ncol_term2;
            Chiso_Jlmp1n_psios_mult = zeros(ncol_term2,S.Nev);
            Chiso_Jlmn_psi_mult = zeros(ncol_term2,S.Nev);

            soindx1 = S.Atom(JJ_a).term2_index_so(1:ncol_term2)+1;
            soindx2 = S.Atom(JJ_a).term2_index_so(1:ncol_term2);

            for img = 1:S.Atom(JJ_a).n_image_rc
                phase_fac = (exp(dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)*fac)));
                Chiso_Jlmp1n_psios_mult = Chiso_Jlmp1n_psios_mult + transpose(bsxfun(@times, conj(S.Atom(JJ_a).rcImage(img).Chiso_mat(:,soindx1)), S.W(S.Atom(JJ_a).rcImage(img).rc_pos))) * S.psi(S.Atom(JJ_a).rcImage(img).rc_pos+S.N,:,kpt) * phase_fac ;
                Chiso_Jlmn_psi_mult = Chiso_Jlmn_psi_mult +         transpose(bsxfun(@times, S.Atom(JJ_a).rcImage(img).Chiso_mat(:,soindx2), S.W(S.Atom(JJ_a).rcImage(img).rc_pos))) * conj(S.psi(S.Atom(JJ_a).rcImage(img).rc_pos,:,kpt)) * conj(phase_fac) ;
            end
            P_nl = P_nl - S.occfac * S.wkpt(kpt) * transpose(S.Atom(JJ_a).term2_gammaso_Jl(1:ncol_term2)) * 2 * real(Chiso_Jlmp1n_psios_mult.*Chiso_Jlmn_psi_mult) * S.occ(:,kpt) ;
        end
        
        % pressure due to scalar relativistic 
        for spinor = 1:S.nspinor
            shift = (spinor-1)*S.N;         % for selecting each spinor, spinor=1 shift = 0, spinor=2, shift = S.N
            nsshift = (spinor-1)*S.tnkpt*(S.spin_typ == 1);

            integral_1 = zeros(S.Atom(JJ_a).angnum,S.Nev);
            integral_2_x = zeros(S.Atom(JJ_a).angnum,S.Nev);
            integral_2_y = zeros(S.Atom(JJ_a).angnum,S.Nev);
            integral_2_z = zeros(S.Atom(JJ_a).angnum,S.Nev);
            
            for img = 1:S.Atom(JJ_a).n_image_rc
                phase_fac = (exp(dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)*fac)));
                ChiW = transpose(bsxfun(@times, conj(S.Atom(JJ_a).rcImage(img).Chi_mat), S.W(S.Atom(JJ_a).rcImage(img).rc_pos)));
                integral_1 = integral_1 + conj(ChiW) * conj(S.psi(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:,kpt)) * conj(phase_fac);
                xr =(S.Atom(JJ_a).rcImage(img).rc_pos_ii-1)*S.dx - S.Atom(JJ_a).rcImage(img).coordinates(1) ;
                yr =(S.Atom(JJ_a).rcImage(img).rc_pos_jj-1)*S.dy - S.Atom(JJ_a).rcImage(img).coordinates(2) ;
                zr =(S.Atom(JJ_a).rcImage(img).rc_pos_kk-1)*S.dz - S.Atom(JJ_a).rcImage(img).coordinates(3) ;
                integral_2_x = integral_2_x + ChiW * ...
                    ((Dpsi_x(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(xr,1,S.Nev)) * phase_fac;
                integral_2_y = integral_2_y + ChiW * ...
                    ((Dpsi_y(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(yr,1,S.Nev)) * phase_fac;
                integral_2_z = integral_2_z + ChiW * ...
                    ((Dpsi_z(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(zr,1,S.Nev)) * phase_fac;
            end

            tf_x = transpose(S.Atom(JJ_a).gamma_Jl) * real(integral_1.*integral_2_x) * S.occ(:,kpt+nsshift);
            tf_y = transpose(S.Atom(JJ_a).gamma_Jl) * real(integral_1.*integral_2_y) * S.occ(:,kpt+nsshift);
            tf_z = transpose(S.Atom(JJ_a).gamma_Jl) * real(integral_1.*integral_2_z) * S.occ(:,kpt+nsshift);
            P_nl = P_nl - 2 * S.occfac * S.wkpt(kpt) * (tf_x + tf_y + tf_z);
        end
        
        % below are all terms related to soc
        if S.Atm(S.Atom(JJ_a).count_typ).pspsoc == 0
            continue;
        end
        
        % pressure due to spin-orbit coupling term 1
        ncol_term1 = S.Atom(JJ_a).ncol_term1;
        soindx = S.Atom(JJ_a).term1_index_so(1:ncol_term1);
        for spinor = 1:S.nspinor
            sigma = (-1)^(spinor-1);
            shift = (spinor-1)*S.N;         % for selecting each spinor, spinor=1 shift = 0, spinor=2, shift = S.N
            integral_1 = zeros(ncol_term1,S.Nev);
            integral_2_x = zeros(ncol_term1,S.Nev);
            integral_2_y = zeros(ncol_term1,S.Nev);
            integral_2_z = zeros(ncol_term1,S.Nev);
            
            for img = 1:S.Atom(JJ_a).n_image_rc
                phase_fac = (exp(dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)*fac)));
                ChisoW = transpose(bsxfun(@times, conj(S.Atom(JJ_a).rcImage(img).Chiso_mat(:,soindx)), S.W(S.Atom(JJ_a).rcImage(img).rc_pos)));
                integral_1 = integral_1 + conj(ChisoW) * conj(S.psi(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:,kpt)) * conj(phase_fac);
                xr =(S.Atom(JJ_a).rcImage(img).rc_pos_ii-1)*S.dx - S.Atom(JJ_a).rcImage(img).coordinates(1) ;
                yr =(S.Atom(JJ_a).rcImage(img).rc_pos_jj-1)*S.dy - S.Atom(JJ_a).rcImage(img).coordinates(2) ;
                zr =(S.Atom(JJ_a).rcImage(img).rc_pos_kk-1)*S.dz - S.Atom(JJ_a).rcImage(img).coordinates(3) ;
                integral_2_x = integral_2_x + ChisoW * ...
                    ((Dpsi_x(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(xr,1,S.Nev)) * phase_fac;
                integral_2_y = integral_2_y + ChisoW * ...
                    ((Dpsi_y(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(yr,1,S.Nev)) * phase_fac;
                integral_2_z = integral_2_z + ChisoW * ...
                    ((Dpsi_z(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(zr,1,S.Nev)) * phase_fac;
            end

            tf_x = transpose(sigma*S.Atom(JJ_a).term1_gammaso_Jl(1:ncol_term1)) * real(integral_1.*integral_2_x) * S.occ(:,kpt);
            tf_y = transpose(sigma*S.Atom(JJ_a).term1_gammaso_Jl(1:ncol_term1)) * real(integral_1.*integral_2_y) * S.occ(:,kpt);
            tf_z = transpose(sigma*S.Atom(JJ_a).term1_gammaso_Jl(1:ncol_term1)) * real(integral_1.*integral_2_z) * S.occ(:,kpt);
            P_nl = P_nl - 2 * S.occfac * S.wkpt(kpt) * (tf_x + tf_y + tf_z);
        end
        
        % pressure due to spin-orbit coupling term 2
        ncol_term2 = S.Atom(JJ_a).ncol_term2;
        for spinor = 1:S.nspinor
            shift = (spinor-1)*S.N;         % for selecting each spinor, spinor=1 shift = 0, spinor=2, shift = S.N
            shift2 = (2-spinor)*S.N;        % for selecting the other spin channel, spinor=1 shift2 = S.N, spinor=2,shift2=0  
            if spinor == 1
                soindx1 = S.Atom(JJ_a).term2_index_so(1:ncol_term2)+1;
                soindx2 = S.Atom(JJ_a).term2_index_so(1:ncol_term2);
            else 
                soindx1 = S.Atom(JJ_a).term2_index_so(1:ncol_term2);
                soindx2 = S.Atom(JJ_a).term2_index_so(1:ncol_term2)+1;
            end
            
            integral_1 = zeros(ncol_term2,S.Nev);
            integral_2_x = zeros(ncol_term2,S.Nev);
            integral_2_y = zeros(ncol_term2,S.Nev);
            integral_2_z = zeros(ncol_term2,S.Nev);
            
            for img = 1:S.Atom(JJ_a).n_image_rc
                phase_fac = (exp(dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)*fac)));
                ChisoW1 = transpose(bsxfun(@times, conj(S.Atom(JJ_a).rcImage(img).Chiso_mat(:,soindx1)), S.W(S.Atom(JJ_a).rcImage(img).rc_pos)));
                ChisoW2 = transpose(bsxfun(@times, S.Atom(JJ_a).rcImage(img).Chiso_mat(:,soindx2), S.W(S.Atom(JJ_a).rcImage(img).rc_pos)));
                
                integral_1 = integral_1 + ChisoW1 * S.psi(S.Atom(JJ_a).rcImage(img).rc_pos+shift2,:,kpt) * phase_fac;

                xr =(S.Atom(JJ_a).rcImage(img).rc_pos_ii-1)*S.dx - S.Atom(JJ_a).rcImage(img).coordinates(1) ;
                yr =(S.Atom(JJ_a).rcImage(img).rc_pos_jj-1)*S.dy - S.Atom(JJ_a).rcImage(img).coordinates(2) ;
                zr =(S.Atom(JJ_a).rcImage(img).rc_pos_kk-1)*S.dz - S.Atom(JJ_a).rcImage(img).coordinates(3) ;
                integral_2_x = integral_2_x + ChisoW2 * ...
                    (conj(Dpsi_x(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(xr,1,S.Nev)) * conj(phase_fac);
                integral_2_y = integral_2_y + ChisoW2 * ...
                    (conj(Dpsi_y(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(yr,1,S.Nev)) * conj(phase_fac);
                integral_2_z = integral_2_z + ChisoW2 * ...
                    (conj(Dpsi_z(S.Atom(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(zr,1,S.Nev)) * conj(phase_fac);
            end

            tf_x = transpose(S.Atom(JJ_a).term2_gammaso_Jl(1:ncol_term2)) * 2*real(integral_1.*integral_2_x) * S.occ(:,kpt);
            tf_y = transpose(S.Atom(JJ_a).term2_gammaso_Jl(1:ncol_term2)) * 2*real(integral_1.*integral_2_y) * S.occ(:,kpt);
            tf_z = transpose(S.Atom(JJ_a).term2_gammaso_Jl(1:ncol_term2)) * 2*real(integral_1.*integral_2_z) * S.occ(:,kpt);
            P_nl = P_nl - S.occfac * S.wkpt(kpt) * (tf_x + tf_y + tf_z);                
        end
    end % end of loop over atoms    
end 

if S.usefock > 1
    pres_exx = evaluateHybridPressure(S);
end

cell_measure = S.Jacb;
if S.BCx == 0
	cell_measure = cell_measure * S.L1;
end
if S.BCy == 0
	cell_measure = cell_measure * S.L2;
end
if S.BCz == 0
	cell_measure = cell_measure * S.L3;
end

if S.usefock > 0
    pressure = -(P_eng + P_elec + P_corr + P_nl + pres_exx + P_nlcc)/(3 * cell_measure);
else
    pressure = -(P_eng + P_elec + P_corr + P_nl + P_nlcc)/(3 * cell_measure);
end

if S.d3Flag == 1 
	pressure = pressure - (S.d3stress(1,1) + S.d3stress(2,2) + S.d3stress(3,3))/3;
end

if (S.vdWDFFlag == 1) || (S.vdWDFFlag == 2)
%     if (S.nspin == 1) % to be removed after finishing svdW-DF
        S = vdWDF_stress(S);
        pressure = pressure - (S.vdWstress(1,1) + S.vdWstress(2,2) + S.vdWstress(3,3))/3;
%     end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    DENSITY MATRIX APPROACH - O(N) PRESSURE CALCULATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function pressure = evaluatePressure(S)

% % Gradient of phi calculation
% Dphi_x = S.grad_1 * S.phi;
% Dphi_y = S.grad_2 * S.phi;
% Dphi_z = S.grad_3 * S.phi;

% % Gradient of rho calculation
% Drho_x = S.grad_1 * S.rho;
% Drho_y = S.grad_2 * S.rho;
% Drho_z = S.grad_3 * S.rho;

% % Component from Energy related terms calculated explicitly

% P_eng = - (sum(S.rho .* (3*S.Vxc + 1.5*S.phi) .* S.W)) + (3*S.Exc) + ...
%         (1.5*sum(S.b .* S.phi .* S.W))  + (3*S.E_corr) + ...
%        - 2*sum(S.W .* S.vdrho_xc .* (S.lapc_T(1,1)*Drho_x.*Drho_x + S.lapc_T(2,2)*Drho_y.*Drho_y + S.lapc_T(3,3)*Drho_z.*Drho_z +...
%         (S.lapc_T(1,2) + S.lapc_T(2,1))*Drho_x.*Drho_y + (S.lapc_T(2,3) + S.lapc_T(3,2))*Drho_y.*Drho_z +...
%         (S.lapc_T(1,3) + S.lapc_T(3,1))*Drho_z.*Drho_x )) + ...
%     (1/4/pi) * sum( (S.lapc_T(1,1)*Dphi_x.*Dphi_x + S.lapc_T(2,2)*Dphi_y.*Dphi_y + S.lapc_T(3,3)*Dphi_z.*Dphi_z +...
%     (S.lapc_T(1,2) + S.lapc_T(2,1))*Dphi_x.*Dphi_y + (S.lapc_T(2,3) + S.lapc_T(3,2))*Dphi_y.*Dphi_z +...
%     (S.lapc_T(1,3) + S.lapc_T(3,1))*Dphi_x.*Dphi_z ) .*S.W ) - 3*S.Eself;

% % Components from gradient terms

% % initialization

% %P_Eself = 3*S.Eself;
% P_elec = 0;
% P_corr = 0;
% dx2 = S.dx*S.dx;
% dy2 = S.dy*S.dy;
% dz2 = S.dz*S.dz;
% dxdy = S.dx*S.dy;
% dydz = S.dy*S.dz;
% dzdx = S.dz*S.dx;
% coeff = S.w2(1) * (S.lapc_T(1,1)/dx2 + S.lapc_T(2,2)/dy2 + S.lapc_T(3,3)/dz2);
% count_typ = 1;
% count_typ_atms = 1;
% for JJ_a = 1:S.n_atm % loop over all the atoms
%     % Load pseudopotential file
%     if count_typ_atms == 1
%         filename = strcat('./Pseudopotentials/Pseudopotential_', S.Atm(count_typ).typ);
%         load(filename)
%     end

%     % Atom position
%     x0 = S.Atoms(JJ_a,1);
%     y0 = S.Atoms(JJ_a,2);
%     z0 = S.Atoms(JJ_a,3);
%     if(S.BC == 2)
%         % Note the S.dx, S.dy, S.dz terms are to ensure the image rb-region overlap w/ fund. domain
%         n_image_xl = floor((S.Atoms(JJ_a,1) + S.Atm(count_typ).rb)/S.L1);
%         n_image_xr = floor((S.L1 - S.Atoms(JJ_a,1)+S.Atm(count_typ).rb-S.dx)/S.L1);
%         n_image_yl = floor((S.Atoms(JJ_a,2) + S.Atm(count_typ).rb)/S.L2);
%         n_image_yr = floor((S.L2 - S.Atoms(JJ_a,2)+S.Atm(count_typ).rb-S.dy)/S.L2);
%         n_image_zl = floor((S.Atoms(JJ_a,3) + S.Atm(count_typ).rb)/S.L3);
%         n_image_zr = floor((S.L3 - S.Atoms(JJ_a,3)+S.Atm(count_typ).rb-S.dz)/S.L3);
%         % Total No. of images of atom JJ_a (including atom JJ_a)
%         n_image_total = (n_image_xl+n_image_xr+1) * (n_image_yl+n_image_yr+1) * (n_image_zl+n_image_zr+1);
%         % Find the coordinates for all the images
%         xx_img = [-n_image_xl : n_image_xr] * S.L1 + x0;
%         yy_img = [-n_image_yl : n_image_yr] * S.L2 + y0;
%         zz_img = [-n_image_zl : n_image_zr] * S.L3 + z0;
%         [XX_IMG_3D,YY_IMG_3D,ZZ_IMG_3D] = meshgrid(xx_img,yy_img,zz_img);
%         XX_IMG_3D = permute(XX_IMG_3D,[2 1 3]);
%         YY_IMG_3D = permute(YY_IMG_3D,[2 1 3]);
%         ZZ_IMG_3D = permute(ZZ_IMG_3D,[2 1 3]);
%     elseif(S.BC == 3)
%         n_image_xl = floor((S.Atoms(JJ_a,1) + S.Atm(count_typ).rb)/S.L1);
%         n_image_xr = floor((S.L1 - S.Atoms(JJ_a,1)+S.Atm(count_typ).rb-S.dx)/S.L1);
%         n_image_yl = floor((S.Atoms(JJ_a,2) + S.Atm(count_typ).rb)/S.L2);
%         n_image_yr = floor((S.L2 - S.Atoms(JJ_a,2)+S.Atm(count_typ).rb-S.dy)/S.L2);
%         % Total No. of images of atom JJ_a (including atom JJ_a)
%         n_image_total = (n_image_xl+n_image_xr+1) * (n_image_yl+n_image_yr+1);
%         % Find the coordinates for all the images
%         xx_img = [-n_image_xl : n_image_xr] * S.L1 + x0;
%         yy_img = [-n_image_yl : n_image_yr] * S.L2 + y0;
%         zz_img = z0;
%         [XX_IMG_3D,YY_IMG_3D,ZZ_IMG_3D] = meshgrid(xx_img,yy_img,zz_img);
%         XX_IMG_3D = permute(XX_IMG_3D,[2 1 3]);
%         YY_IMG_3D = permute(YY_IMG_3D,[2 1 3]);
%         ZZ_IMG_3D = permute(ZZ_IMG_3D,[2 1 3]);
%     elseif(S.BC == 4)
%         n_image_zl = floor((S.Atoms(JJ_a,3) + S.Atm(count_typ).rb)/S.L3);
%         n_image_zr = floor((S.L3 - S.Atoms(JJ_a,3)+S.Atm(count_typ).rb-S.dz)/S.L3);

%         % Total No. of images of atom JJ_a (including atom JJ_a)
%         n_image_total = n_image_zl+n_image_zr+1;
%         % Find the coordinates for all the images
%         xx_img = x0;
%         yy_img = y0;
%         zz_img = [-n_image_zl : n_image_zr] * S.L3 + z0;
%         [XX_IMG_3D,YY_IMG_3D,ZZ_IMG_3D] = meshgrid(xx_img,yy_img,zz_img);
%         XX_IMG_3D = permute(XX_IMG_3D,[2 1 3]);
%         YY_IMG_3D = permute(YY_IMG_3D,[2 1 3]);
%         ZZ_IMG_3D = permute(ZZ_IMG_3D,[2 1 3]);
%     end

%     % Loop over all image(s) of atom JJ_a (including atom JJ_a)
%     for count_image = 1:n_image_total

%         % Atom position of the image
%         x0_i = XX_IMG_3D(count_image);
%         y0_i = YY_IMG_3D(count_image);
%         z0_i = ZZ_IMG_3D(count_image);

%         % Indices of closest grid point to atom
%         pos_ii = round(x0_i / S.dx) + 1;
%         pos_jj = round(y0_i / S.dy) + 1;
%         pos_kk = round(z0_i / S.dz) + 1;

%         % Starting and ending indices of b-region
%         ii_s = pos_ii - ceil(S.Atm(count_typ).rb/S.dx+0.5);
%         ii_e = pos_ii + ceil(S.Atm(count_typ).rb/S.dx+0.5);
%         jj_s = pos_jj - ceil(S.Atm(count_typ).rb/S.dy+0.5);
%         jj_e = pos_jj + ceil(S.Atm(count_typ).rb/S.dy+0.5);
%         kk_s = pos_kk - ceil(S.Atm(count_typ).rb/S.dz+0.5);
%         kk_e = pos_kk + ceil(S.Atm(count_typ).rb/S.dz+0.5);

%         if(S.BC == 2)
%             % For periodic systems, find the overlap of rb-region and the fund. domain
%             ii_s = max(ii_s,1);
%             ii_e = min(ii_e,S.Nx);
%             jj_s = max(jj_s,1);
%             jj_e = min(jj_e,S.Ny);
%             kk_s = max(kk_s,1);
%             kk_e = min(kk_e,S.Nz);
%         elseif(S.BC == 3)
%             ii_s = max(ii_s,1);
%             ii_e = min(ii_e,S.Nx);
%             jj_s = max(jj_s,1);
%             jj_e = min(jj_e,S.Ny);
%             isInside = (kk_s>1) && (kk_e<S.Nz);
%             assert(isInside,'ERROR: Atom too close to boundary in z-direction for b calculation');
%         elseif(S.BC == 4)
%             isInside = (ii_s>1) && (ii_e<S.Nx) && (jj_s>1) && (jj_e<S.Ny);
%             assert(isInside,'ERROR: Atom too close to boundary for b calculation');
%             kk_s = max(kk_s,1);
%             kk_e = min(kk_e,S.Nz);
%         end

%         xx = (ii_s-2*S.FDn-1:ii_e+2*S.FDn-1)*S.dx - x0_i;
%         yy = (jj_s-2*S.FDn-1:jj_e+2*S.FDn-1)*S.dy - y0_i;
%         zz = (kk_s-2*S.FDn-1:kk_e+2*S.FDn-1)*S.dz - z0_i;
%         [XX_3D,YY_3D,ZZ_3D] = meshgrid(xx,yy,zz);
%         XX_3D = permute(XX_3D,[2 1 3]);
%         YY_3D = permute(YY_3D,[2 1 3]);
%         ZZ_3D = permute(ZZ_3D,[2 1 3]);
%         dd = sqrt(S.metric_T(1,1)*XX_3D.^2 + S.metric_T(1,2)*(XX_3D.*YY_3D) + S.metric_T(1,3)*(XX_3D.*ZZ_3D) + ...
%             S.metric_T(2,1)*(YY_3D.*XX_3D) + S.metric_T(2,2)*YY_3D.^2 + S.metric_T(2,3)*(YY_3D.*ZZ_3D) + ...
%             S.metric_T(3,1)*(ZZ_3D.*XX_3D) + S.metric_T(3,2)*(ZZ_3D.*YY_3D) + S.metric_T(3,3)*ZZ_3D.^2) ;

%         % Pseudopotential at grid points through interpolation
%         V_PS = interp1(S.Atm(count_typ).r_grid_vloc, S.Atm(count_typ).r_grid_vloc.*S.Atm(count_typ).Vloc, dd, 'spline');
%         V_PS = V_PS./dd;
%         V_PS(dd<S.Atm(count_typ).r_grid_vloc(2)) = S.Atm(count_typ).Vloc(1);
%         % Reference potential at grid points
%         rc_ref = S.rc_ref; % WARNING: Might need smaller if pseudocharges overlap
%         V_PS_ref = zeros(size(dd));
%         I_ref = dd<rc_ref;
%         V_PS_ref(~I_ref) = -(S.Atm(count_typ).Z)./dd(~I_ref);
%         V_PS_ref(I_ref) = -S.Atm(count_typ).Z*(9*dd(I_ref).^7-30*rc_ref*dd(I_ref).^6 ...
%             +28*rc_ref*rc_ref*dd(I_ref).^5-14*(rc_ref^5)*dd(I_ref).^2+12*rc_ref^7)/(5*rc_ref^8);

%         % Pseudocharge density
%         bJ = zeros(size(V_PS));
%         bJ_ref = zeros(size(V_PS));
%         II = 1+S.FDn : size(V_PS,1)-S.FDn;
%         JJ = 1+S.FDn : size(V_PS,2)-S.FDn;
%         KK = 1+S.FDn : size(V_PS,3)-S.FDn;
%         bJ(II,JJ,KK) = coeff * V_PS(II,JJ,KK);
%         bJ_ref(II,JJ,KK) = coeff * V_PS_ref(II,JJ,KK);
%         for p = 1:S.FDn
%             bJ(II,JJ,KK) = bJ(II,JJ,KK) + S.w2(p+1)*S.lapc_T(1,1)/dx2 * (V_PS(II+p,JJ,KK) + V_PS(II-p,JJ,KK)) + ...
%                 S.w2(p+1)*S.lapc_T(2,2)/dy2 * (V_PS(II,JJ+p,KK) + V_PS(II,JJ-p,KK)) + ...
%                 S.w2(p+1)*S.lapc_T(3,3)/dz2 * (V_PS(II,JJ,KK+p) + V_PS(II,JJ,KK-p));
%             bJ_ref(II,JJ,KK) = bJ_ref(II,JJ,KK) + S.w2(p+1)*S.lapc_T(1,1)/dx2 * (V_PS_ref(II+p,JJ,KK) + V_PS_ref(II-p,JJ,KK)) + ...
%                 S.w2(p+1)*S.lapc_T(2,2)/dy2 * (V_PS_ref(II,JJ+p,KK) + V_PS_ref(II,JJ-p,KK)) + ...
%                 S.w2(p+1)*S.lapc_T(3,3)/dz2 * (V_PS_ref(II,JJ,KK+p) + V_PS_ref(II,JJ,KK-p));
%             for q = 1:S.FDn
%                 bJ(II,JJ,KK) = bJ(II,JJ,KK) + S.w1(p+1)*S.w1(q+1)*(S.lapc_T(1,2) + S.lapc_T(2,1))/dxdy * ( V_PS(II+q,JJ+p,KK) - ...
%                     V_PS(II-q,JJ+p,KK) - V_PS(II+q,JJ-p,KK) + V_PS(II-q,JJ-p,KK) ) + ...
%                     S.w1(p+1)*S.w1(q+1)*(S.lapc_T(2,3) + S.lapc_T(3,2))/dydz * ( V_PS(II,JJ+q,KK+p) - ...
%                     V_PS(II,JJ-q,KK+p) - V_PS(II,JJ+q,KK-p) + V_PS(II,JJ-q,KK-p) ) + ...
%                     S.w1(p+1)*S.w1(q+1)*(S.lapc_T(1,3) + S.lapc_T(3,1))/dzdx * ( V_PS(II+q,JJ,KK+p) - ...
%                     V_PS(II-q,JJ,KK+p) - V_PS(II+q,JJ,KK-p) + V_PS(II-q,JJ,KK-p) ) ;
%                 bJ_ref(II,JJ,KK) = bJ_ref(II,JJ,KK) + S.w1(p+1)*S.w1(q+1)*(S.lapc_T(1,2) + S.lapc_T(2,1))/dxdy * ( V_PS_ref(II+q,JJ+p,KK) - ...
%                     V_PS_ref(II-q,JJ+p,KK) - V_PS_ref(II+q,JJ-p,KK) + V_PS_ref(II-q,JJ-p,KK) ) + ...
%                     S.w1(p+1)*S.w1(q+1)*(S.lapc_T(2,3) + S.lapc_T(3,2))/dydz * ( V_PS_ref(II,JJ+q,KK+p) - ...
%                     V_PS_ref(II,JJ-q,KK+p) - V_PS_ref(II,JJ+q,KK-p) + V_PS_ref(II,JJ-q,KK-p) ) + ...
%                     S.w1(p+1)*S.w1(q+1)*(S.lapc_T(1,3) + S.lapc_T(3,1))/dzdx * ( V_PS_ref(II+q,JJ,KK+p) - ...
%                     V_PS_ref(II-q,JJ,KK+p) - V_PS_ref(II+q,JJ,KK-p) + V_PS_ref(II-q,JJ,KK-p) ) ;
%             end
%         end
%         bJ = (-1/(4*pi))*bJ;
%         bJ_ref = (-1/(4*pi))*bJ_ref;

%         % Calculate the gradient of pseudocharges
%         dbJ_x = zeros(size(V_PS)); dbJ_y = zeros(size(V_PS)); dbJ_z = zeros(size(V_PS));
%         dbJ_ref_x = zeros(size(V_PS)); dbJ_ref_y = zeros(size(V_PS)); dbJ_ref_z = zeros(size(V_PS));
%         dVJ_x = zeros(size(V_PS)); dVJ_y = zeros(size(V_PS)); dVJ_z = zeros(size(V_PS));
%         dVJ_ref_x = zeros(size(V_PS)); dVJ_ref_y = zeros(size(V_PS)); dVJ_ref_z = zeros(size(V_PS));

%         II = 1+2*S.FDn : size(V_PS,1)-2*S.FDn;
%         JJ = 1+2*S.FDn : size(V_PS,2)-2*S.FDn;
%         KK = 1+2*S.FDn : size(V_PS,3)-2*S.FDn;

%         for p = 1:S.FDn
%             dbJ_x(II,JJ,KK) = dbJ_x(II,JJ,KK) + S.w1(p+1)/S.dx*(bJ(II+p,JJ,KK)-bJ(II-p,JJ,KK));
%             dbJ_y(II,JJ,KK) = dbJ_y(II,JJ,KK) + S.w1(p+1)/S.dy*(bJ(II,JJ+p,KK)-bJ(II,JJ-p,KK));
%             dbJ_z(II,JJ,KK) = dbJ_z(II,JJ,KK) + S.w1(p+1)/S.dz*(bJ(II,JJ,KK+p)-bJ(II,JJ,KK-p));
%             dbJ_ref_x(II,JJ,KK) = dbJ_ref_x(II,JJ,KK) + S.w1(p+1)/S.dx*(bJ_ref(II+p,JJ,KK)-bJ_ref(II-p,JJ,KK));
%             dbJ_ref_y(II,JJ,KK) = dbJ_ref_y(II,JJ,KK) + S.w1(p+1)/S.dy*(bJ_ref(II,JJ+p,KK)-bJ_ref(II,JJ-p,KK));
%             dbJ_ref_z(II,JJ,KK) = dbJ_ref_z(II,JJ,KK) + S.w1(p+1)/S.dz*(bJ_ref(II,JJ,KK+p)-bJ_ref(II,JJ,KK-p));
%             dVJ_x(II,JJ,KK) = dVJ_x(II,JJ,KK) + S.w1(p+1)/S.dx*(V_PS(II+p,JJ,KK)-V_PS(II-p,JJ,KK));
%             dVJ_y(II,JJ,KK) = dVJ_y(II,JJ,KK) + S.w1(p+1)/S.dy*(V_PS(II,JJ+p,KK)-V_PS(II,JJ-p,KK));
%             dVJ_z(II,JJ,KK) = dVJ_z(II,JJ,KK) + S.w1(p+1)/S.dz*(V_PS(II,JJ,KK+p)-V_PS(II,JJ,KK-p));
%             dVJ_ref_x(II,JJ,KK) = dVJ_ref_x(II,JJ,KK) + S.w1(p+1)/S.dx*(V_PS_ref(II+p,JJ,KK)-V_PS_ref(II-p,JJ,KK));
%             dVJ_ref_y(II,JJ,KK) = dVJ_ref_y(II,JJ,KK) + S.w1(p+1)/S.dy*(V_PS_ref(II,JJ+p,KK)-V_PS_ref(II,JJ-p,KK));
%             dVJ_ref_z(II,JJ,KK) = dVJ_ref_z(II,JJ,KK) + S.w1(p+1)/S.dz*(V_PS_ref(II,JJ,KK+p)-V_PS_ref(II,JJ,KK-p));

%         end

%         % Calculate local force and correction force components
%         [II_rb,JJ_rb,KK_rb] = meshgrid(ii_s:ii_e,jj_s:jj_e,kk_s:kk_e);
%         II_rb = permute(II_rb,[2,1,3]);
%         JJ_rb = permute(JJ_rb,[2,1,3]);
%         KK_rb = permute(KK_rb,[2,1,3]);
%         Rowcount_rb = (KK_rb-1)*S.Nx*S.Ny + (JJ_rb-1)*S.Nx + II_rb;
%         [xr,yr,zr] = meshgrid((ii_s-1:ii_e-1)*S.dx - x0_i,(jj_s-1:jj_e-1)*S.dy - y0_i,(kk_s-1:kk_e-1)*S.dz - z0_i) ;
%         xr = permute(xr,[2 1 3]);
%         yr = permute(yr,[2 1 3]);
%         zr = permute(zr,[2 1 3]);

%         P_elec = P_elec + sum(sum(sum( dVJ_x(II,JJ,KK) .* xr .* ( - 0.5 * bJ(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
%         P_elec = P_elec + sum(sum(sum( dVJ_y(II,JJ,KK) .* yr .* ( - 0.5 * bJ(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
%         P_elec = P_elec + sum(sum(sum( dVJ_z(II,JJ,KK) .* zr .* ( - 0.5 * bJ(II,JJ,KK) ) .* S.W(Rowcount_rb) )));

%         P_elec = P_elec + sum(sum(sum( dbJ_x(II,JJ,KK) .* xr .* ( S.phi(Rowcount_rb) -0.5 * V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
%         P_elec = P_elec + sum(sum(sum( dbJ_y(II,JJ,KK) .* yr .* ( S.phi(Rowcount_rb) -0.5 * V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));
%         P_elec = P_elec + sum(sum(sum( dbJ_z(II,JJ,KK) .* zr .* ( S.phi(Rowcount_rb) -0.5 * V_PS(II,JJ,KK) ) .* S.W(Rowcount_rb) )));

%         P_corr = P_corr + 0.5 * sum(sum(sum( ( dbJ_x(II,JJ,KK) .* ( S.V_c(Rowcount_rb) + V_PS(II,JJ,KK) ) + ...
%             dbJ_ref_x(II,JJ,KK) .* ( S.V_c(Rowcount_rb) - V_PS_ref(II,JJ,KK) ) + ...
%             dVJ_ref_x(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ_ref(II,JJ,KK)) - ...
%             dVJ_x(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ(II,JJ,KK)) ) .* xr .* S.W(Rowcount_rb) )));

%         P_corr = P_corr + 0.5 * sum(sum(sum( ( dbJ_y(II,JJ,KK) .* ( S.V_c(Rowcount_rb) + V_PS(II,JJ,KK) ) + ...
%             dbJ_ref_y(II,JJ,KK) .* ( S.V_c(Rowcount_rb) - V_PS_ref(II,JJ,KK) ) + ...
%             dVJ_ref_y(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ_ref(II,JJ,KK)) - ...
%             dVJ_y(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ(II,JJ,KK)) ) .* yr .* S.W(Rowcount_rb) )));

%         P_corr = P_corr + 0.5 * sum(sum(sum( ( dbJ_z(II,JJ,KK) .* ( S.V_c(Rowcount_rb) + V_PS(II,JJ,KK) ) + ...
%             dbJ_ref_z(II,JJ,KK) .* ( S.V_c(Rowcount_rb) - V_PS_ref(II,JJ,KK) ) + ...
%             dVJ_ref_z(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ_ref(II,JJ,KK)) - ...
%             dVJ_z(II,JJ,KK) .* ( S.b(Rowcount_rb) + S.b_ref(Rowcount_rb) - bJ(II,JJ,KK)) ) .* zr .* S.W(Rowcount_rb) )));

%     end


%     % Check if same type of atoms are over
%     if count_typ_atms == S.Atm(count_typ).n_atm_typ
%         count_typ_atms = 1;
%         count_typ = count_typ + 1;
%     else
%         count_typ_atms = count_typ_atms + 1;
%     end

% end % end of loop over atoms


% %**********************************************************************
% %*                   Calculate Kinetic + nonlocal Pressure                  *
% %**********************************************************************


% P_nl = 0;
% fac = [S.L1 S.L2 S.L3];
% Pmat = zeros(S.Nx*S.Ny*S.Nz,S.Nx*S.Ny*S.Nz,S.tnkpt);
% parfor kpt = 1:S.tnkpt
%     for i = 1: S.Nev
%         Pmat(:,:,kpt) = Pmat(:,:,kpt) + S.psi(:,i,kpt)*transpose(conj(S.psi(:,i,kpt)))*S.occ(i,kpt);
%     end
% end

% parfor kpt = 1:S.tnkpt
%     kpt_vec = S.kptgrid(kpt,:);

%     % Density matrix formation
% %    for i = 1: S.Nev
% %        Pmat(:,:,kpt) = Pmat(:,:,kpt) + S.psi(:,i,kpt)*S.psi(:,i,kpt)'*S.occ(i,kpt);
% %    end

%     Grad_1 = S.grad_T(1,1)*blochGradient(S,kpt_vec,1) + S.grad_T(2,1)*blochGradient(S,kpt_vec,2) + S.grad_T(3,1)*blochGradient(S,kpt_vec,3);
%     Grad_2 = S.grad_T(1,2)*blochGradient(S,kpt_vec,1) + S.grad_T(2,2)*blochGradient(S,kpt_vec,2) + S.grad_T(3,2)*blochGradient(S,kpt_vec,3);
%     Grad_3 = S.grad_T(1,3)*blochGradient(S,kpt_vec,1) + S.grad_T(2,3)*blochGradient(S,kpt_vec,2) + S.grad_T(3,3)*blochGradient(S,kpt_vec,3);

%     % Gradient of density matrix
%     GPmat_x = blochGradient(S,kpt_vec,1) * Pmat(:,:,kpt);
%     GPmat_y = blochGradient(S,kpt_vec,2) * Pmat(:,:,kpt);
%     GPmat_z = blochGradient(S,kpt_vec,3) * Pmat(:,:,kpt);

%     % Laplacian of the density matrix
%     LPmat_1 = Grad_1*Pmat(:,:,kpt)*Grad_1;
%     LPmat_2 = Grad_2*Pmat(:,:,kpt)*Grad_2;
%     LPmat_3 = Grad_3*Pmat(:,:,kpt)*Grad_3;

%     % Kinetic pressure
%     P_eng = P_eng + real(2*S.wkpt(kpt)*trace((LPmat_1 + LPmat_2 + LPmat_3).*diag(S.W)));

%     % Non-local pressure
%     count_typ = 1;
%     count_typ_atms = 1;
%     for JJ_a = 1:S.n_atm % loop over all atoms

%         for l = 0:S.Atom(count_typ).lmax
%             if l == S.Atom(count_typ).lloc
%                 continue;
%             end
%             for m = -l:l
%                 Chi_X_mult1 = zeros(1,S.Nx*S.Ny*S.Nz);

%                 for img = 1:S.Atom(JJ_a).n_image_rc
%                     phase_fac = (exp(1i*2*pi*dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)./fac)));
%                     Chi_X_mult1 = Chi_X_mult1 + transpose(conj((S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat(:,m+l+1) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)))) * Pmat(S.Atom(JJ_a).rcImage(img).rc_pos,:,kpt) * phase_fac ;
%                 end

%                 for img = 1:S.Atom(JJ_a).n_image_rc
%                     phase_fac = (exp(-1i*2*pi*dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)./fac)));
%                     P_nl = P_nl - 6 * S.wkpt(kpt) * S.Atom(count_typ).gamma_Jl(l+1) * Chi_X_mult1(1,S.Atom(JJ_a).rcImage(img).rc_pos) * ...
%                         (S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat(:,m+l+1).* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)) * phase_fac ;
%                 end

%                 Chi_X_mult1 = 0*Chi_X_mult1;

%                 for img = 1:S.Atom(JJ_a).n_image_rc
%                     phase_fac = (exp(1i*2*pi*dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)./fac)));
%                     xr =(S.Atom(JJ_a).rcImage(img).rc_pos_ii-1)*S.dx - S.Atom(JJ_a).rcImage(img).coordinates(1) ;
%                     yr =(S.Atom(JJ_a).rcImage(img).rc_pos_jj-1)*S.dy - S.Atom(JJ_a).rcImage(img).coordinates(2) ;
%                     zr =(S.Atom(JJ_a).rcImage(img).rc_pos_kk-1)*S.dz - S.Atom(JJ_a).rcImage(img).coordinates(3) ;
%                     Chi_X_mult1 = Chi_X_mult1 + (transpose(conj(S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat(:,m+l+1)) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos).*xr) * ...
%                         GPmat_x(S.Atom(JJ_a).rcImage(img).rc_pos,:) + transpose(conj(S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat(:,m+l+1)) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos).*yr) * ...
%                         GPmat_y(S.Atom(JJ_a).rcImage(img).rc_pos,:) + transpose(conj(S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat(:,m+l+1)) .* S.W(S.Atom(JJ_a).rcImage(img).rc_pos).*zr) * ...
%                         GPmat_z(S.Atom(JJ_a).rcImage(img).rc_pos,:)) * phase_fac;
%                 end

%                 for img = 1:S.Atom(JJ_a).n_image_rc
%                     phase_fac = (exp(-1i*2*pi*dot(kpt_vec,(S.Atoms(JJ_a,:)-S.Atom(JJ_a).rcImage(img).coordinates)./fac)));
%                     P_nl = P_nl - 4 * S.wkpt(kpt) * S.Atom(count_typ).gamma_Jl(l+1) * real((Chi_X_mult1(1,S.Atom(JJ_a).rcImage(img).rc_pos)) * ...
%                         (S.Atom(JJ_a).Chi(l+1).rcImage(img).Chi_mat(:,m+l+1).* S.W(S.Atom(JJ_a).rcImage(img).rc_pos)) * phase_fac) ;
%                 end
%             end
%         end

%         % Check if same type of atoms are over
%         if count_typ_atms == S.Atm(count_typ).n_atm_typ
%             count_typ_atms = 1;
%             count_typ = count_typ + 1;
%         else
%             count_typ_atms = count_typ_atms + 1;
%         end

%     end % end of loop over atoms
% end
% pressure = -(P_eng + P_elec + P_corr + P_nl)/(3 * cell_measure);
