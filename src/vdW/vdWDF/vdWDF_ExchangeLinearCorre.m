function [S, ecPW, v_cPW] = vdWDF_ExchangeLinearCorre(S, XC)
% @file    vdWDF_ExchangeLinearCorre.m
% @brief   This file contains the functions for computing exchange and linear 
%          correlation energy density and potential of vdW-DF1 and vdW-DF2.
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
    % exchange part
    if S.NLCC_flag 
        rho = S.rho+S.rho_Tilde_at;
    else 
        rho = S.rho;
    end
	rho(rho < S.xc_rhotol) = S.xc_rhotol;
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

    if S.vdWDFFlag == 1 % vdWDF1: Zhang-Yang revPBE
        XC.kappa_pbe = 1.245; % Zhang-Yang revPBE
        XC.kappa = XC.kappa_pbe;
        XC.mu_divkappa = XC.mu/XC.kappa_pbe;
%         % Arrays
% 	    rho_updn = rho/2.0;
% 	    rho_updnm1_3 = rho_updn.^(-XC.third);
% 	    rhom1_3 = XC.twom1_3 * rho_updnm1_3;

% 	    rhotot_inv = rhom1_3 .* rhom1_3 .* rhom1_3;
% % 	    rhotmo6 = sqrt(rhom1_3);
% % 	    rhoto6 = rho .* rhom1_3 .* rhom1_3 .* rhotmo6;

% 	    rhomot = rho_updnm1_3;
% 	    ex_lsd = -XC.threefourth_divpi * XC.sixpi2_1_3 * (rhomot .* rhomot .* rho_updn);

% 	    %        Perdew-Burke-Ernzerhof GGA, exchange part, changed kappa

% 	    rho_inv = rhomot .* rhomot .* rhomot;
% 	    coeffss = (1.0/4.0) * XC.sixpi2m1_3 * XC.sixpi2m1_3 * (rho_inv .* rho_inv .* rhomot .* rhomot);
% 	    ss = (sigma/4.0) .* coeffss;
% 	    divss = 1.0./(1.0 + XC.mu_divkappa * ss);
% 	    dfxdss = XC.mu * (divss .* divss);
% 	    %d2fxdss2 = -XC.mu * 2.0 * XC.mu_divkappa * (divss .* divss .* divss);
% 	    fx = 1.0 + XC.kappa * (1.0 - divss);
% 	    ex_gga = ex_lsd .* fx;
% 	    dssdn = (-8.0/3.0) * (ss .* rho_inv);
% 	    dfxdn = dfxdss .* dssdn;
% 	    v_x = ex_lsd .* ((4.0/3.0) * fx + rho_updn .* dfxdn);

% 	    dssdg = 2.0 * coeffss;
% 	    dfxdg = dfxdss .* dssdg;
% 	    dvxdgrho1 = 0.5*ex_lsd .* rho_updn .* dfxdg;
% 	    ex = ex_gga .* rho_updn;

% 	    %        If non spin-polarized, treat spin down contribution now, similar to spin up
% 	    ex = ex * 2.0;
%         ex = ex .* rhotot_inv;

		% below is QuantumEspresso PBE code
        agrho = sigma.^0.5;
        kf = (3*pi^2)^(1/3) * rho.^(1/3);
        dsg = 0.5 ./ kf;
        s1 = agrho .* dsg ./ rho;
        s2 = s1 .* s1;
        f1 = s2 * XC.mu / XC.kappa;
        f2 = 1.0 + f1;
        f3 = XC.kappa ./ f2;
        fx = 1 + XC.kappa - f3; % here is different from QE, 1 added

        exunif = - 3/(4*pi) * kf;
        ex = exunif .* fx;

        dxunif = exunif * 1/3;
        ds = - 4/3 * s1;

        dfx1 = f2 .* f2;
        dfx = 2.0 * XC.mu * s1 ./ dfx1;

        v_x= ex + dxunif .* fx + exunif .* dfx .* ds;
        dvxdgrho1 = exunif .* dfx .* dsg ./ agrho;
    elseif S.vdWDFFlag == 2 % vdWDF2: GGA revised PW86
        a = 1.851;
        b = 17.33;
        c = 0.163;
        s_prefactor = 6.18733545256027; % 2*(3\pi^2)^(1/3)
        Ax = -0.738558766382022; % -3/4 * (3/pi)^(1/3)
        four_thirds = 4.0/3.0;
        
        grad_rho = sigma.^0.5;
        s = grad_rho ./ (s_prefactor*rho.^four_thirds);
        s_2 = s.*s;
        s_3 = s_2.*s;
        s_4 = s_3.*s;
        s_5 = s_3.*s_2;
        s_6 = s_5.*s;
        fs = (1.0 + a*s_2 + b*s_4 + c*s_6).^(1.0/15.0);
        ex = Ax * rho.^(1.0/3.0) .* fs; % \epsilon_x, not n\epsilon_x
        df_ds = (1.0./(15.0*fs.^14.0)) .* (2.0*a*s + 4.0*b*s_3 + 6.0*c*s_5);
        v_x = Ax*four_thirds * (rho.^(1.0/3.0) .*fs - grad_rho./(s_prefactor*rho).*df_ds);
        dvxdgrho1 = Ax * df_ds./(s_prefactor*grad_rho);
    end
    % correlation part, LDA PW91
	p = 1 ;
	A = 0.031091 ;
	alpha1 = 0.21370 ;
	beta1 = 7.5957 ;
	beta2 = 3.5876 ;
	beta3 = 1.6382 ;
	beta4 = 0.49294 ;

	% correlation potential
    rs = (0.75./(pi*rho)).^(1/3) ;
	G1 = log(1+1./(2*A*(beta1*(rs.^0.5) + beta2*rs + beta3*(rs.^1.5) + beta4*(rs.^(p+1)))));
    ec = -2.0 * A * (1 + alpha1 * rs) .* G1;
	v_c = (-2*A*(1+alpha1*rs)).*G1 ...
		- (rs/3).*(-2*A*alpha1*G1 ...
		- ((-2*A*(1+alpha1*rs)).*(A*( beta1*(rs.^-0.5)+ 2*beta2 + 3*beta3*(rs.^0.5) + 2*(p+1)*beta4*(rs.^p) ))) ...
		./((2*A*( beta1*(rs.^0.5) + beta2*rs + beta3*(rs.^1.5) + beta4*(rs.^(p+1)) ) ) ...
		.*(2*A*( beta1*(rs.^0.5) + beta2*rs + beta3*(rs.^1.5) + beta4*(rs.^(p+1)) ) )+(2*A*( beta1*(rs.^0.5) + beta2*rs + beta3*(rs.^1.5) + beta4*(rs.^(p+1)) ) )) ) ;
    % summary
    S.e_xc = ex + ec;
    ecPW = ec;
    v_cPW = v_c;
	S.dvxcdgrho = dvxdgrho1;

    if S.cell_typ ~= 2
		S.Vxc = v_x + v_c - S.grad_1 * (S.dvxcdgrho.*drho_1) - S.grad_2 * (S.dvxcdgrho.*drho_2) - S.grad_3 * (S.dvxcdgrho.*drho_3);
	else
		S.Vxc = v_x + v_c - ( S.lapc_T(1,1)*S.grad_1*(S.dvxcdgrho.*drho_1) + S.lapc_T(2,2)*S.grad_2*(S.dvxcdgrho.*drho_2) + S.lapc_T(3,3)*S.grad_3*(S.dvxcdgrho.*drho_3) +...
						 S.lapc_T(2,1)*S.grad_1*(S.dvxcdgrho.*drho_2) + S.lapc_T(2,1)*S.grad_2*(S.dvxcdgrho.*drho_1) + S.lapc_T(3,2)*S.grad_2*(S.dvxcdgrho.*drho_3) +...
						 S.lapc_T(3,2)*S.grad_3*(S.dvxcdgrho.*drho_2) + S.lapc_T(3,1)*S.grad_1*(S.dvxcdgrho.*drho_3) + S.lapc_T(3,1)*S.grad_3*(S.dvxcdgrho.*drho_1) );
    end
end