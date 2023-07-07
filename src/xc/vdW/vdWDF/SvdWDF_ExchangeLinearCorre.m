function [S, ecPW, v_cPW] = SvdWDF_ExchangeLinearCorre(S, XC)
% @file    SvdWDF_ExchangeLinearCorre.m
% @brief   This file contains the functions for computing exchange and linear 
%          correlation energy density and potential of spin-polarized vdW-DF1 and vdW-DF2.
% @authors Boqin Zhang <bzhang376@gatech.edu>
%          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% Reference:
% Dion, Max, Henrik Rydberg, Elsebeth Schröder, David C. Langreth, and Bengt I. Lundqvist. 
% "Van der Waals density functional for general geometries." 
% Physical review letters 92, no. 24 (2004): 246401.
% Román-Pérez, Guillermo, and José M. Soler. 
% "Efficient implementation of a van der Waals density functional: application to double-wall carbon nanotubes." 
% Physical review letters 103, no. 9 (2009): 096102.
% Lee, Kyuho, Éamonn D. Murray, Lingzhu Kong, Bengt I. Lundqvist, and David C. Langreth. 
% "Higher-accuracy van der Waals density functional." Physical Review B 82, no. 8 (2010): 081101.
% Thonhauser, T., S. Zuluaga, C. A. Arter, K. Berland, E. Schröder, and P. Hyldgaard. 
% "Spin signature of nonlocal correlation binding in metal-organic frameworks." 
% Physical review letters 115, no. 13 (2015): 136402.
% Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
% ==============================================================================================
    rho = S.rho;
    if S.NLCC_flag 
        rho(:,2) = rho(:,2)+S.rho_Tilde_at * 0.5;
        rho(:,3) = rho(:,3)+S.rho_Tilde_at * 0.5;
    end
    rho(rho < S.xc_rhotol) = S.xc_rhotol;
	rho(:,1) = rho(:,2) + rho(:,3);
    
    drho_1 = S.grad_1 * rho;
	drho_2 = S.grad_2 * rho;
	drho_3 = S.grad_3 * rho;

    if S.cell_typ ~= 2
		sigma = drho_1.*drho_1 + drho_2.*drho_2 + drho_3.*drho_3;
    else
		sigma = (S.lapc_T(1,1)*drho_1.*drho_1 + S.lapc_T(2,2)*drho_2.*drho_2 + S.lapc_T(3,3)*drho_3.*drho_3 +...
				 S.lapc_T(1,2)*drho_1.*drho_2 + S.lapc_T(2,3)*drho_2.*drho_3 + S.lapc_T(1,3)*drho_3.*drho_1 ) ; % grad_rho . grad_rho
    end
	sigma(sigma < 1E-14) = 1E-14;
    
    % Arrays
    rho_updnm1_3 = rho(:,2:3).^(-XC.third);
    rhom1_3 = rho(:,1).^(-XC.third);
    rhotot_inv = rhom1_3.^3;
    zeta = (rho(:,2) - rho(:,3)) .* rhotot_inv;
    zetp = 1.0 + zeta * XC.alpha_zeta;
    zetm = 1.0 - zeta * XC.alpha_zeta;
    zetpm1_3 = zetp.^(-XC.third);       
    zetmm1_3 = zetm.^(-XC.third);

    rhotmo6 = sqrt(rhom1_3);
    rhoto6 = rho(:,1) .* rhom1_3 .* rhom1_3 .* rhotmo6;
    
    % -----------------------------------------------------------------------
	% First take care of the exchange part of the functional
    if S.vdWDFFlag == 1 % vdWDF1: Zhang-Yang revPBE
        XC.kappa_pbe = 1.245; % Zhang-Yang revPBE
        XC.kappa = XC.kappa_pbe;
        XC.mu_divkappa = XC.mu/XC.kappa_pbe;

        rhomot = rho_updnm1_3;
        ex_lsd = -XC.threefourth_divpi * XC.sixpi2_1_3 * (rhomot .* rhomot .* rho(:,2:3));

        rho_inv = rhomot .* rhomot .* rhomot;
        coeffss = (1.0/4.0) * XC.sixpi2m1_3 * XC.sixpi2m1_3 * (rho_inv .* rho_inv .* rhomot .* rhomot);
        ss = sigma(:,2:3) .* coeffss;
    
        divss = 1.0./(1.0 + XC.mu_divkappa * ss);
        dfxdss = XC.mu * (divss .* divss);
    
        fx = 1.0 + XC.kappa * (1.0 - divss);
        ex_gga = ex_lsd .* fx;
        dssdn = (-8.0/3.0) * (ss .* rho_inv);
        dfxdn = dfxdss .* dssdn;
        v_x = ex_lsd .* ((4.0/3.0) * fx + rho(:,2:3) .* dfxdn);

        dssdg = 2.0 * coeffss;
        dfxdg = dfxdss .* dssdg; 
        dvxdgrho1 = ex_lsd .* rho(:,2:3) .* dfxdg;
        ex = sum(ex_gga .* rho(:,2:3),2) .* rhotot_inv;
    elseif S.vdWDFFlag == 2 % vdWDF2: GGA revised PW86
        a = 1.851;
        b = 17.33;
        c = 0.163;
        s_prefactor = 6.18733545256027; % 2*(3\pi^2)^(1/3)
        Ax = -0.738558766382022; % -3/4 * (3/pi)^(1/3)
        four_thirds = 4.0/3.0;
        
        grad_rho = sigma(:,2:3).^0.5;
        s = grad_rho ./ (2^(1/3) * s_prefactor*rho(:, 2:3).^four_thirds);
        s_2 = s.*s;
        s_3 = s_2.*s;
        s_4 = s_3.*s;
        s_5 = s_3.*s_2;
        s_6 = s_5.*s;
        fs = (1.0 + a*s_2 + b*s_4 + c*s_6).^(1.0/15.0);
        ex = Ax*2^(1/3) * sum(rho(:, 2:3).^(1.0/3.0 + 1) .* fs, 2) ./ rho(:, 1); % \epsilon_x, not n\epsilon_x
        df_ds = (1.0./(15.0*fs.^14.0)) .* (2.0*a*s + 4.0*b*s_3 + 6.0*c*s_5);
        v_x = Ax*four_thirds * (2^(1/3)*rho(:, 2:3).^(1.0/3.0) .*fs - grad_rho./(s_prefactor*rho(:, 2:3)).*df_ds);
        dvxdgrho1 = Ax * df_ds./(s_prefactor*grad_rho);
    end
    
    % -----------------------------------------------------------------------------
	% Then takes care of the LSD correlation part of the functional

	rs = XC.rsfac * rhom1_3;
	sqr_rs = XC.sq_rsfac * rhotmo6;
	rsm1_2 = XC.sq_rsfac_inv * rhoto6;

	% Formulas A6-A8 of PW92LSD

	ec0_q0 = -2.0 * XC.ec0_aa * (1.0 + XC.ec0_a1 * rs);
	ec0_q1 = 2.0 * XC.ec0_aa *(XC.ec0_b1 * sqr_rs + XC.ec0_b2 * rs + XC.ec0_b3 * rs .* sqr_rs + XC.ec0_b4 * rs .* rs);
	ec0_q1p = XC.ec0_aa * (XC.ec0_b1 * rsm1_2 + 2.0 * XC.ec0_b2 + 3.0 * XC.ec0_b3 * sqr_rs + 4.0 * XC.ec0_b4 * rs);
	ec0_den = 1.0./(ec0_q1 .* ec0_q1 + ec0_q1);
	ec0_log = -log(ec0_q1 .* ec0_q1 .* ec0_den);
	ecrs0 = ec0_q0 .* ec0_log;
	decrs0_drs = -2.0 * XC.ec0_aa * XC.ec0_a1 * ec0_log - ec0_q0 .* ec0_q1p .* ec0_den;

	mac_q0 = -2.0 * XC.mac_aa * (1.0 + XC.mac_a1 * rs);
	mac_q1 = 2.0 * XC.mac_aa * (XC.mac_b1 * sqr_rs + XC.mac_b2 * rs + XC.mac_b3 * rs .* sqr_rs + XC.mac_b4 * rs .* rs);
	mac_q1p = XC.mac_aa * (XC.mac_b1 * rsm1_2 + 2 * XC.mac_b2 + 3 * XC.mac_b3 * sqr_rs + 4 * XC.mac_b4 * rs);
	mac_den = 1.0./(mac_q1 .* mac_q1 + mac_q1);
	mac_log = -log( mac_q1 .* mac_q1 .* mac_den );
	macrs = mac_q0 .* mac_log;
	dmacrs_drs = -2.0 * XC.mac_aa * XC.mac_a1 * mac_log - mac_q0 .* mac_q1p .* mac_den;

	%zeta = (rho(:,2) - rho(:,3)) .* rhotot_inv;
	ec1_q0 = -2.0 * XC.ec1_aa * (1.0 + XC.ec1_a1 * rs);
	ec1_q1 = 2.0 * XC.ec1_aa * (XC.ec1_b1 * sqr_rs + XC.ec1_b2 * rs + XC.ec1_b3 * rs .* sqr_rs + XC.ec1_b4 * rs .* rs);
	ec1_q1p = XC.ec1_aa * (XC.ec1_b1 * rsm1_2 + 2 * XC.ec1_b2 + 3 * XC.ec1_b3 * sqr_rs + 4 * XC.ec1_b4 * rs);
	ec1_den = 1.0./(ec1_q1 .* ec1_q1 + ec1_q1);
	ec1_log = -log( ec1_q1 .* ec1_q1 .* ec1_den );
	ecrs1 = ec1_q0 .* ec1_log;
	decrs1_drs = -2.0 * XC.ec1_aa * XC.ec1_a1 * ec1_log - ec1_q0 .* ec1_q1p .* ec1_den;
	
	% XC.alpha_zeta is introduced in order to remove singularities for fully polarized systems.
	zetp_1_3 = (1.0 + zeta * XC.alpha_zeta) .* (zetpm1_3.^2);
	zetm_1_3 = (1.0 - zeta * XC.alpha_zeta) .* (zetmm1_3.^2);

	f_zeta = ( (1.0 + zeta * XC.alpha_zeta2) .* zetp_1_3 + (1.0 - zeta * XC.alpha_zeta2) .* zetm_1_3 - 2.0 ) * XC.factf_zeta;
	fp_zeta = ( zetp_1_3 - zetm_1_3 ) * XC.factfp_zeta;
	zeta4 = zeta.^4;

	gcrs = ecrs1 - ecrs0 + macrs * XC.fsec_inv;
	ec = ecrs0 + f_zeta .* (zeta4 .* gcrs - macrs * XC.fsec_inv);
	dgcrs_drs = decrs1_drs - decrs0_drs + dmacrs_drs * XC.fsec_inv;
	decrs_drs = decrs0_drs + f_zeta .* (zeta4 .* dgcrs_drs - dmacrs_drs * XC.fsec_inv);
	dfzeta4_dzeta = 4.0 * zeta.^3 .* f_zeta + fp_zeta .* zeta4;
	decrs_dzeta = dfzeta4_dzeta .* gcrs - fp_zeta .* macrs * XC.fsec_inv;

	S.e_xc = ex + ec;
	vxcadd = ec - rs * XC.third .* decrs_drs - zeta .* decrs_dzeta;
    ecPW = ec;
    v_cPW = [vxcadd + decrs_dzeta, vxcadd - decrs_dzeta];
	v_xc(:,1) = v_x(:,1) + vxcadd + decrs_dzeta;
	v_xc(:,2) = v_x(:,2) + vxcadd - decrs_dzeta;
	
    S.dvxcdgrho = zeros(S.N,3);
	S.dvxcdgrho(:, 2:3) = dvxdgrho1; % 1st column should be Vc for grad rho; in vdWDF it is zero vector
    if S.cell_typ ~= 2
		Vxc_temp = S.grad_1 * (dvxdgrho1.*drho_1(:, 2:3)) + S.grad_2 * (dvxdgrho1.*drho_2(:, 2:3)) + S.grad_3 * (dvxdgrho1.*drho_3(:, 2:3));
	else
		Vxc_temp =  S.lapc_T(1,1)*S.grad_1*(dvxdgrho1.*drho_1(:, 2:3)) + S.lapc_T(2,2)*S.grad_2*(dvxdgrho1.*drho_2(:, 2:3)) + S.lapc_T(3,3)*S.grad_3*(dvxdgrho1.*drho_3(:, 2:3)) +...
                    S.lapc_T(2,1)*S.grad_1*(dvxdgrho1.*drho_2(:, 2:3)) + S.lapc_T(2,1)*S.grad_2*(dvxdgrho1.*drho_1(:, 2:3)) + S.lapc_T(3,2)*S.grad_2*(dvxdgrho1.*drho_3(:, 2:3)) +...
					S.lapc_T(3,2)*S.grad_3*(dvxdgrho1.*drho_2(:, 2:3)) + S.lapc_T(3,1)*S.grad_1*(dvxdgrho1.*drho_3(:, 2:3)) + S.lapc_T(3,1)*S.grad_3*(dvxdgrho1.*drho_1(:, 2:3));
    end
    S.Vxc = v_xc - Vxc_temp;
end