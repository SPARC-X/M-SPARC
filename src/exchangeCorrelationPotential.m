function [S] = exchangeCorrelationPotential(S)
% @ brief    Calculates exchange correlation potential(LDA, GGA) and 
%             exchange correlation energy density(GGA). Purdew-Wang LDA 
%             and PBE GGA are implemented for both spin polarized and
%             non-spin polarized calculations.
% @ authors
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Qimen Xu <qimenxu@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%==========================================================================

%--------------------------------------------------------------------------    
% Parameters
%
XC = struct;
XC.alpha_zeta2 = 1.0 - 1.0e-6; XC.alpha_zeta = 1.0 - 1.0e-6; % ABINIT
%XC.alpha_zeta2 = 1.0; XC.alpha_zeta = 1.0; %LIBXC
if (strcmp(S.XC, 'GGA_PBEsol'))
    XC.beta = 0.046;
    XC.mu = 10.0/81.0;
else
    XC.beta = 0.066725;
    %XC.beta = 0.06672455060314922;
    XC.mu = 0.2195149727645171;
end
XC.fsec_inv = 1.0/1.709921;
XC.kappa_pbe = 0.804;
XC.rsfac = 0.6203504908994000;
XC.kappa = XC.kappa_pbe;
XC.mu_divkappa_pbe = XC.mu/XC.kappa_pbe;
XC.mu_divkappa = XC.mu_divkappa_pbe;

XC.ec0_aa = 0.031091; XC.ec1_aa = 0.015545; XC.mac_aa = 0.016887; % ABINIT
%XC.ec0_aa = 0.0310907; XC.ec1_aa = 0.01554535; XC.mac_aa = 0.0168869; % LIBXC
XC.ec0_a1 = 0.21370;  XC.ec1_a1 = 0.20548;  XC.mac_a1 = 0.11125;
XC.ec0_b1 = 7.5957;  XC.ec1_b1 = 14.1189;  XC.mac_b1 = 10.357;
XC.ec0_b2 = 3.5876;   XC.ec1_b2 = 6.1977;   XC.mac_b2 = 3.6231;
XC.ec0_b3 = 1.6382;   XC.ec1_b3 = 3.3662;   XC.mac_b3 = 0.88026;
XC.ec0_b4 = 0.49294;  XC.ec1_b4 = 0.62517;  XC.mac_b4 = 0.49671;

% Constants
XC.piinv = 1.0/pi;
XC.third = 1.0/3.0;
XC.twom1_3 = 2.0^(-XC.third);
XC.sixpi2_1_3 = (6.0 * pi^2)^XC.third;
XC.sixpi2m1_3 = 1.0/XC.sixpi2_1_3;
XC.threefourth_divpi = (3.0/4.0) * XC.piinv;
XC.gamma = (1.0 - log(2.0)) * XC.piinv^2;
XC.gamma_inv = 1.0/XC.gamma;
%beta_gamma = XC.beta * XC.gamma_inv;
XC.factf_zeta = 1.0/(2.0^(4.0/3.0) - 2.0);
XC.factfp_zeta = 4.0/3.0 * XC.factf_zeta * XC.alpha_zeta2;
XC.coeff_tt = 1.0/(4.0 * 4.0 * XC.piinv * (3.0 * pi^2)^XC.third);
XC.sq_rsfac = sqrt(XC.rsfac);
XC.sq_rsfac_inv = 1.0/XC.sq_rsfac;


if S.nspin == 1
	if S.xc == 0
		S = LDA_PW(S);
	elseif S.xc == 1
		S = LDA_PZ(S);
	elseif S.xc == 2
		S = GGA_PBE(S,XC);
	elseif (S.xc == -102) || (S.xc == -108)
        S = Calculate_Vxc_vdWExchangeLinearCorre(S, XC);
        S = Calculate_nonLinearCorr_E_V_vdWDF(S);
	end
else
	if S.xc == 0
		S = LSDA_PW(S,XC);
	elseif S.xc == 1
		S = LSDA_PZ(S,XC);
	elseif S.xc == 2
		S = GSGA_PBE(S,XC);
	end
end

end
%--------------------------------------------------------------------------

function [S] = LDA_PW(S)
    if S.NLCC_flag 
        rho = S.rho+S.rho_Tilde_at;
    else 
        rho = S.rho;
    end
	rho(rho < S.xc_rhotol) = S.xc_rhotol;
	% correlation parameters
	p = 1 ;
	A = 0.031091 ;
	alpha1 = 0.21370 ;
	beta1 = 7.5957 ;
	beta2 = 3.5876 ;
	beta3 = 1.6382 ;
	beta4 = 0.49294 ;

	% exchange parameters
	C3 = 0.9847450218427;

	% exchange-correlation potential
	%rho = rho+(1e-50) ;
	S.Vxc = (0.75./(pi*rho)).^(1/3) ;
	S.Vxc = - C3*(rho.^(1/3)) + (-2*A*(1+alpha1*S.Vxc)).*log(1+1./(2*A*(beta1*(S.Vxc.^0.5) + beta2*S.Vxc + beta3*(S.Vxc.^1.5) + beta4*(S.Vxc.^(p+1))))) ...
		- (S.Vxc/3).*(-2*A*alpha1*log(1+1./(2*A*( beta1*(S.Vxc.^0.5) + beta2*S.Vxc + beta3*(S.Vxc.^1.5) + beta4*(S.Vxc.^(p+1))))) ...
		- ((-2*A*(1+alpha1*S.Vxc)).*(A*( beta1*(S.Vxc.^-0.5)+ 2*beta2 + 3*beta3*(S.Vxc.^0.5) + 2*(p+1)*beta4*(S.Vxc.^p) ))) ...
		./((2*A*( beta1*(S.Vxc.^0.5) + beta2*S.Vxc + beta3*(S.Vxc.^1.5) + beta4*(S.Vxc.^(p+1)) ) ) ...
		.*(2*A*( beta1*(S.Vxc.^0.5) + beta2*S.Vxc + beta3*(S.Vxc.^1.5) + beta4*(S.Vxc.^(p+1)) ) )+(2*A*( beta1*(S.Vxc.^0.5) + beta2*S.Vxc + beta3*(S.Vxc.^1.5) + beta4*(S.Vxc.^(p+1)) ) )) ) ;
	S.dvxcdgrho = zeros(size(S.Vxc));
end
%--------------------------------------------------------------------------

function [S] = LDA_PZ(S) 
    if S.NLCC_flag 
        rho = S.rho+S.rho_Tilde_at;
    else 
        rho = S.rho;
    end
	rho(rho < S.xc_rhotol) = S.xc_rhotol;
	% correlation parameters
	A = 0.0311;
	B = -0.048;
	C = 0.002;
	D = -0.0116;
	gamma1 = -0.1423;
	beta1 = 1.0529;
	beta2 = 0.3334;
	% exchange parameters
	C3 = 0.9847450218427;
	% exchange-correlation potential
	%rho = rho+(1e-50);
	S.Vxc = (0.75./(pi*rho)).^(1/3);
	islessthan1 = (S.Vxc < 1.0);
	S.Vxc(islessthan1) = log(S.Vxc(islessthan1)).*(A+(2.0/3.0)*C*S.Vxc(islessthan1)) ...
		+ (B-(1.0/3.0)*A) + (1.0/3.0)*(2.0*D-C)* S.Vxc(islessthan1);
	S.Vxc(~islessthan1) = (gamma1 + (7.0/6.0)*gamma1*beta1*sqrt(S.Vxc(~islessthan1)) ...
		+ (4.0/3.0)*gamma1*beta2*S.Vxc(~islessthan1))./(1+beta1*sqrt(S.Vxc(~islessthan1))+beta2*S.Vxc(~islessthan1)).^2;
	S.Vxc = S.Vxc - C3*(rho.^(1/3));
	isRhoZero = (abs(rho) < 1e-15);
	S.Vxc(isRhoZero) = 0;
	%rho = rho-(1e-50) ;
	S.dvxcdgrho = zeros(size(rho,1),1);
end
%------------------------------------------------------------------------------------------------------------------------------    

function [S] = GGA_PBE(S,XC)
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
	% -------------------------------------------------------------------------------------------------------------------
	%                                                     Calling inbuilt Libxc functions
	% -------------------------------------------------------------------------------------------------------------------

%     func_id_x = 101; % Perdew-Burke-Ernzerhof GGA functional with exchange only
%     %mex -I/home/users/asharma424/libxc/src -I/home/users/asharma424/libxc -L/opt/etsf/lib -lxc -lm /home/users/asharma424/libxc/examples/xc.c ;
%     fprintf('\n\n');
%     mex -I/home/users/asharma424/Downloads/Softwares/libxc-4.3.3/src -I/home/users/asharma424/Downloads/Softwares/libxc-4.3.3 -L/usr/lib64 -lxc -lm /home/users/asharma424/Downloads/Softwares/libxc-4.3.3/examples/xc.c ;
%     [e_x, v_x, dvxdgrho] = xc(func_id_x, rho, sigma);% e_x = exchange energy per unit particle, v_x = d(rho*e_x)/drho = e_x + rho*de_x/drho, dvxdgrho = d(rho*e_x)/dsigma = rho*de_x/dsigma + 0
%     
%     func_id_c = 130; % Perdew-Burke-Ernzerhof GGA functional with correlation only
%     mex -I/home/users/asharma424/Downloads/Softwares/libxc-4.3.3/src -I/home/users/asharma424/Downloads/Softwares/libxc-4.3.3 -L/usr/lib64 -lxc -lm /home/users/asharma424/Downloads/Softwares/libxc-4.3.3/examples/xc.c ;
%     [e_c, v_c, dvcdgrho] = xc(func_id_c, rho, sigma);% e_c = correlation energy per unit particle, v_c = d(rho*e_c)/drho = e_c + rho*de_c/drho, dvcdgrho = d(rho*e_c)/dsigma = rho*de_c/dsigma + 0
%     fprintf('\n\n');
%     S.e_xc = e_x + e_c ;
%     v_xc = v_x + v_c ;
%     S.dvxcdgrho = dvxdgrho + dvcdgrho ;
% % 
%     S.Vxc = v_xc - 2*( S.lapc_T(1,1)*grad_1*(S.dvxcdgrho.*drho_1) + S.lapc_T(2,2)*grad_2*(S.dvxcdgrho.*drho_2) + S.lapc_T(3,3)*grad_3*(S.dvxcdgrho.*drho_3) +...
%                      S.lapc_T(2,1)*grad_1*(S.dvxcdgrho.*drho_2) + S.lapc_T(2,1)*grad_2*(S.dvxcdgrho.*drho_1) + S.lapc_T(3,2)*grad_2*(S.dvxcdgrho.*drho_3) +...
%                      S.lapc_T(3,2)*grad_3*(S.dvxcdgrho.*drho_2) + S.lapc_T(3,1)*grad_1*(S.dvxcdgrho.*drho_3) + S.lapc_T(3,1)*grad_3*(S.dvxcdgrho.*drho_1) );     
	% --------------------------------------------------------------------------------------------------------------------------
	%                                         Direct computation
	%                                         (taken from ABINIT)
	% --------------------------------------------------------------------------------------------------------------------------
   
	
	% Arrays
	rho_updn = rho/2.0;
	rho_updnm1_3 = rho_updn.^(-XC.third);
	rhom1_3 = XC.twom1_3 * rho_updnm1_3;

	rhotot_inv = rhom1_3 .* rhom1_3 .* rhom1_3;
	rhotmo6 = sqrt(rhom1_3);
	rhoto6 = rho .* rhom1_3 .* rhom1_3 .* rhotmo6;

	%        -----------------------------------------------------------------------
	%        First take care of the exchange part of the functional

	rhomot = rho_updnm1_3;
	ex_lsd = -XC.threefourth_divpi * XC.sixpi2_1_3 * (rhomot .* rhomot .* rho_updn);

	%        Perdew-Burke-Ernzerhof GGA, exchange part

	rho_inv = rhomot .* rhomot .* rhomot;
	coeffss = (1.0/4.0) * XC.sixpi2m1_3 * XC.sixpi2m1_3 * (rho_inv .* rho_inv .* rhomot .* rhomot);
	ss = (sigma/4.0) .* coeffss;
    
    if (strcmp(S.XC,'GGA_PBE') || strcmp(S.XC,'GGA_PBEsol'))
        divss = 1.0./(1.0 + XC.mu_divkappa * ss);
        dfxdss = XC.mu * (divss .* divss);
    elseif (strcmp(S.XC,'GGA_RPBE'))
        divss = exp(-XC.mu_divkappa * ss);
        dfxdss = XC.mu * divss;
    end
	
	fx = 1.0 + XC.kappa * (1.0 - divss);
	ex_gga = ex_lsd .* fx;
	dssdn = (-8.0/3.0) * (ss .* rho_inv);
	dfxdn = dfxdss .* dssdn;
	v_xc = ex_lsd .* ((4.0/3.0) * fx + rho_updn .* dfxdn);

	dssdg = 2.0 * coeffss;
	dfxdg = dfxdss .* dssdg;
	dvxcdgrho1 = ex_lsd .* rho_updn .* dfxdg;
	exc = ex_gga .* rho_updn;

	%        If non spin-polarized, treat spin down contribution now, similar to spin up

	exc = exc * 2.0;
	S.e_xc = exc .* rhotot_inv;

	%        -----------------------------------------------------------------------------
	%        Then takes care of the LSD correlation part of the functional

	rs = XC.rsfac * rhom1_3;
	sqr_rs = XC.sq_rsfac * rhotmo6;
	rsm1_2 = XC.sq_rsfac_inv * rhoto6;

	%        Formulas A6-A8 of PW92LSD

	ec0_q0 = -2.0 * XC.ec0_aa * (1.0 + XC.ec0_a1 * rs);
	ec0_q1 = 2.0 * XC.ec0_aa *(XC.ec0_b1 * sqr_rs + XC.ec0_b2 * rs + XC.ec0_b3 * rs .* sqr_rs + XC.ec0_b4 * rs .* rs);
	ec0_q1p = XC.ec0_aa * (XC.ec0_b1 * rsm1_2 + 2.0 * XC.ec0_b2 + 3.0 * XC.ec0_b3 * sqr_rs + 4.0 * XC.ec0_b4 * rs);
	ec0_den = 1.0./(ec0_q1 .* ec0_q1 + ec0_q1);
	ec0_log = -log(ec0_q1 .* ec0_q1 .* ec0_den);
	ecrs0 = ec0_q0 .* ec0_log;
	decrs0_drs = -2.0 * XC.ec0_aa * XC.ec0_a1 * ec0_log - ec0_q0 .* ec0_q1p .* ec0_den;

	ecrs = ecrs0;
	decrs_drs = decrs0_drs;
	%decrs_dzeta = 0.0;
	%zeta = 0.0;

	%        Add LSD correlation functional to GGA exchange functional
	S.e_xc = S.e_xc + ecrs;
	v_xc = v_xc + ecrs - (rs/3.0) .* decrs_drs;

	%        -----------------------------------------------------------------------------
	%        Eventually add the GGA correlation part of the PBE functional
	%        Note : the computation of the potential in the spin-unpolarized
	%        case could be optimized much further. Other optimizations are left to do.

	%phi_zeta = 1.0;
	%phip_zeta = 0.0;
	phi_zeta_inv = 1.0;
	%phi_logder = 0.0;
	phi3_zeta = 1.0;
	gamphi3inv = XC.gamma_inv;
	%phipp_zeta = (-2.0/9.0) * XC.alpha_zeta * XC.alpha_zeta;

	%        From ec to bb
	bb = ecrs * gamphi3inv;
	dbb_drs = decrs_drs * gamphi3inv;
	%dbb_dzeta = gamphi3inv * (decrs_dzeta - 3.0 * ecrs * phi_logder);

	%        From bb to cc
	exp_pbe = exp(-bb);
	cc = 1.0./(exp_pbe - 1.0);
	dcc_dbb = cc .* cc .* exp_pbe;
	dcc_drs = dcc_dbb .* dbb_drs;
	%dcc_dzeta = dcc_dbb .* dbb_dzeta;

	%        From cc to aa
	coeff_aa = XC.beta * XC.gamma_inv * phi_zeta_inv * phi_zeta_inv;
	aa = coeff_aa * cc;
	daa_drs = coeff_aa * dcc_drs;
	%daa_dzeta = -2.0 * aa * phi_logder + coeff_aa * dcc_dzeta;

	%        Introduce tt : do not assume that the spin-dependent gradients are collinear
	grrho2 = sigma;
	dtt_dg = 2.0 * rhotot_inv .* rhotot_inv .* rhom1_3 * XC.coeff_tt;
	%        Note that tt is (the t variable of PBE divided by phi) squared
	tt = 0.5 * grrho2 .* dtt_dg;

	%        Get xx from aa and tt
	xx = aa .* tt;
	dxx_drs = daa_drs .* tt;
	%dxx_dzeta = daa_dzeta .* tt;
	dxx_dtt = aa;

	%        From xx to pade
	pade_den = 1.0./(1.0 + xx .* (1.0 + xx));
	pade = (1.0 + xx) .* pade_den;
	dpade_dxx = -xx .* (2.0 + xx) .* (pade_den.^2);
	dpade_drs = dpade_dxx .* dxx_drs;
	dpade_dtt = dpade_dxx .* dxx_dtt;
	%dpade_dzeta = dpade_dxx .* dxx_dzeta;

	%        From pade to qq
	coeff_qq = tt * phi_zeta_inv * phi_zeta_inv;
	qq = coeff_qq .* pade;
	dqq_drs = coeff_qq .* dpade_drs;
	dqq_dtt = pade * phi_zeta_inv * phi_zeta_inv + coeff_qq .* dpade_dtt;
	%dqq_dzeta = coeff_qq .* (dpade_dzeta - 2.0 * pade * phi_logder);

	%        From qq to rr
	arg_rr = 1.0 + XC.beta * XC.gamma_inv * qq;
	div_rr = 1.0./arg_rr;
	rr = XC.gamma * log(arg_rr);
	drr_dqq = XC.beta * div_rr;
	drr_drs = drr_dqq .* dqq_drs;
	drr_dtt = drr_dqq .* dqq_dtt;
	%drr_dzeta = drr_dqq .* dqq_dzeta;

	%        From rr to hh
	hh = phi3_zeta * rr;
	dhh_drs = phi3_zeta * drr_drs;
	dhh_dtt = phi3_zeta * drr_dtt;
	%dhh_dzeta = phi3_zeta * (drr_dzeta + 3.0 * rr * phi_logder);

	%        The GGA correlation energy is added
	S.e_xc = S.e_xc + hh;

	%        From hh to the derivative of the energy wrt the density
	drhohh_drho = hh - XC.third * rs .* dhh_drs - (7.0/3.0) * tt .* dhh_dtt; %- zeta * dhh_dzeta 
	v_xc = v_xc + drhohh_drho;

	%        From hh to the derivative of the energy wrt to the gradient of the
	%        density, divided by the gradient of the density
	%        (The v3.3 definition includes the division by the norm of the gradient)

	S.dvxcdgrho = 0.5 * dvxcdgrho1 + (rho .* dtt_dg .* dhh_dtt);

	if S.cell_typ ~= 2
		S.Vxc = v_xc - S.grad_1 * (S.dvxcdgrho.*drho_1) - S.grad_2 * (S.dvxcdgrho.*drho_2) - S.grad_3 * (S.dvxcdgrho.*drho_3);
	else
		S.Vxc = v_xc - ( S.lapc_T(1,1)*S.grad_1*(S.dvxcdgrho.*drho_1) + S.lapc_T(2,2)*S.grad_2*(S.dvxcdgrho.*drho_2) + S.lapc_T(3,3)*S.grad_3*(S.dvxcdgrho.*drho_3) +...
						 S.lapc_T(2,1)*S.grad_1*(S.dvxcdgrho.*drho_2) + S.lapc_T(2,1)*S.grad_2*(S.dvxcdgrho.*drho_1) + S.lapc_T(3,2)*S.grad_2*(S.dvxcdgrho.*drho_3) +...
						 S.lapc_T(3,2)*S.grad_3*(S.dvxcdgrho.*drho_2) + S.lapc_T(3,1)*S.grad_1*(S.dvxcdgrho.*drho_3) + S.lapc_T(3,1)*S.grad_3*(S.dvxcdgrho.*drho_1) );
	end
end
%--------------------------------------------------------------------------


function [S] = Calculate_Vxc_vdWExchangeLinearCorre(S, XC)
    % exchange part
    rho = S.rho;
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
	S.dvxcdgrho = dvxdgrho1;

    if S.cell_typ ~= 2
		S.Vxc = v_x + v_c - S.grad_1 * (S.dvxcdgrho.*drho_1) - S.grad_2 * (S.dvxcdgrho.*drho_2) - S.grad_3 * (S.dvxcdgrho.*drho_3);
	else
		S.Vxc = v_x + v_c - ( S.lapc_T(1,1)*S.grad_1*(S.dvxcdgrho.*drho_1) + S.lapc_T(2,2)*S.grad_2*(S.dvxcdgrho.*drho_2) + S.lapc_T(3,3)*S.grad_3*(S.dvxcdgrho.*drho_3) +...
						 S.lapc_T(2,1)*S.grad_1*(S.dvxcdgrho.*drho_2) + S.lapc_T(2,1)*S.grad_2*(S.dvxcdgrho.*drho_1) + S.lapc_T(3,2)*S.grad_2*(S.dvxcdgrho.*drho_3) +...
						 S.lapc_T(3,2)*S.grad_3*(S.dvxcdgrho.*drho_2) + S.lapc_T(3,1)*S.grad_1*(S.dvxcdgrho.*drho_3) + S.lapc_T(3,1)*S.grad_3*(S.dvxcdgrho.*drho_1) );
    end
end

% This function is for calculating non-local correlation part of van der Waals functional.
% Reference: 
% M. Dion, H. Rydberg, E. Schroeder, D.C. Langreth, and B.I. Lundqvist, Phys. Rev. Lett. 92, 246401 (2004)
% G. Roman-Perez and J.M. Soler, Phys. Rev. Lett. 103, 096102 (2009)
function [S] = Calculate_nonLinearCorr_E_V_vdWDF(S)
    S = vdWDF_getQ0onGrid(S);
    S = vdWDF_splineInterpolation_energy(S); % compute the vector u for potential and vdW energy
    S = vdWDF_uGenerate_Potential(S);
    S.Vxc = S.Vxc + S.vdWpotential;
end
%--------------------------------------------------------------------------


function [S] = LSDA_PW(S,XC)
	rho = S.rho;
    if S.NLCC_flag 
        rho(:,2) = rho(:,2)+S.rho_Tilde_at * 0.5;
        rho(:,3) = rho(:,3)+S.rho_Tilde_at * 0.5;
    end
	rho(rho < S.xc_rhotol) = S.xc_rhotol;
	rho(:,1) = rho(:,2) + rho(:,3);
	% Arrays
	%rho = rho+(1e-50) ;
	rho_updnm1_3 = rho(:,2:3).^(-XC.third);
	rhom1_3 = rho(:,1).^(-XC.third);
	rhotot_inv = rhom1_3.^3;
	zeta = (rho(:,2) - rho(:,3)) .* rhotot_inv; % Check whether it is rho_up-rho_dn or rho_dn-rho_up
	zetp = 1 + zeta * XC.alpha_zeta;
	zetm = 1 - zeta * XC.alpha_zeta;
	zetpm1_3 = zetp.^(-XC.third);
	zetmm1_3 = zetm.^(-XC.third);

	rhotmo6 = sqrt(rhom1_3);
	rhoto6 = rho(:,1) .* rhom1_3 .* rhom1_3 .* rhotmo6;

	% -----------------------------------------------------------------------
	% First take care of the exchange part of the functional

	rhomot = rho_updnm1_3;
	ex_lsd = -XC.threefourth_divpi * XC.sixpi2_1_3 * (rhomot .* rhomot .* rho(:,2:3));

	v_xc = (4/3) * ex_lsd;
	exc = sum(ex_lsd .* rho(:,2:3),2);
	S.e_xc = exc .* rhotot_inv;

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
	ecrs = ecrs0 + f_zeta .* (zeta4 .* gcrs - macrs * XC.fsec_inv);
	dgcrs_drs = decrs1_drs - decrs0_drs + dmacrs_drs * XC.fsec_inv;
	decrs_drs = decrs0_drs + f_zeta .* (zeta4 .* dgcrs_drs - dmacrs_drs * XC.fsec_inv);
	dfzeta4_dzeta = 4.0 * zeta.^3 .* f_zeta + fp_zeta .* zeta4;
	decrs_dzeta = dfzeta4_dzeta .* gcrs - fp_zeta .* macrs * XC.fsec_inv;

	S.e_xc = S.e_xc + ecrs;
	vxcadd = ecrs - rs * XC.third .* decrs_drs - zeta .* decrs_dzeta;
	v_xc(:,1) = v_xc(:,1) + vxcadd + decrs_dzeta;
	v_xc(:,2) = v_xc(:,2) + vxcadd - decrs_dzeta;
	S.Vxc = v_xc;
	
	S.dvxcdgrho = zeros(S.N,3);
end
%--------------------------------------------------------------------------


function [S] = LSDA_PZ(S,XC)
	error('LSDA_PZ is not implemented!');
end


%--------------------------------------------------------------------------
function [S] = GSGA_PBE(S,XC)
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

	rhomot = rho_updnm1_3;
	ex_lsd = -XC.threefourth_divpi * XC.sixpi2_1_3 * (rhomot .* rhomot .* rho(:,2:3));

	rho_inv = rhomot .* rhomot .* rhomot;
	coeffss = (1.0/4.0) * XC.sixpi2m1_3 * XC.sixpi2m1_3 * (rho_inv .* rho_inv .* rhomot .* rhomot);
	ss = sigma(:,2:3) .* coeffss;
    
    if (strcmp(S.XC,'GGA_PBE') || strcmp(S.XC,'GGA_PBEsol'))
        divss = 1.0./(1.0 + XC.mu_divkappa * ss);
        dfxdss = XC.mu * (divss .* divss);
    elseif (strcmp(S.XC,'GGA_RPBE'))
        divss = exp(-XC.mu_divkappa * ss);
        dfxdss = XC.mu * divss;
    end
    
	fx = 1.0 + XC.kappa * (1.0 - divss);
	ex_gga = ex_lsd .* fx;
	dssdn = (-8.0/3.0) * (ss .* rho_inv);
	dfxdn = dfxdss .* dssdn;
	v_xc = ex_lsd .* ((4.0/3.0) * fx + rho(:,2:3) .* dfxdn);

	dssdg = 2.0 * coeffss;
	dfxdg = dfxdss .* dssdg; 
	dvxcdgrho1 = ex_lsd .* rho(:,2:3) .* dfxdg;
	exc = sum(ex_gga .* rho(:,2:3),2);

	S.e_xc = exc .* rhotot_inv;

	% -----------------------------------------------------------------------------
	% Then takes care of the LSD correlation part of the functional

	rs = XC.rsfac * rhom1_3;
	sqr_rs = XC.sq_rsfac * rhotmo6;
	rsm1_2 = XC.sq_rsfac_inv * rhoto6;

	%        Formulas A6-A8 of PW92LSD

	ec0_q0 = -2.0 * XC.ec0_aa * (1.0 + XC.ec0_a1 * rs);
	ec0_q1 = 2.0 * XC.ec0_aa *(XC.ec0_b1 * sqr_rs + XC.ec0_b2 * rs + XC.ec0_b3 * rs .* sqr_rs + XC.ec0_b4 * rs .* rs);
	ec0_q1p = XC.ec0_aa * (XC.ec0_b1 * rsm1_2 + 2.0 * XC.ec0_b2 + 3.0 * XC.ec0_b3 * sqr_rs + 4.0 * XC.ec0_b4 * rs);
	ec0_den = 1.0./(ec0_q1 .* ec0_q1 + ec0_q1);
	ec0_log = -log(ec0_q1 .* ec0_q1 .* ec0_den);
	%ec0_log = log( 1.0 + 1.0./ec0_q1);
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
	ecrs = ecrs0 + f_zeta .* (zeta4 .* gcrs - macrs * XC.fsec_inv);
	
	%ecrs = ecrs0 + f_zeta .* (-macrs.* (1.0-zeta4) * XC.fsec_inv + (ecrs1-ecrs0) .* zeta4);
	dgcrs_drs = decrs1_drs - decrs0_drs + dmacrs_drs * XC.fsec_inv;
	decrs_drs = decrs0_drs + f_zeta .* (zeta4 .* dgcrs_drs - dmacrs_drs * XC.fsec_inv);
	%decrs_drs = decrs0_drs + f_zeta .* (-dmacrs_drs .* (1.0 - zeta4) * XC.fsec_inv + (decrs1_drs-decrs0_drs) .* zeta4);
	dfzeta4_dzeta = 4.0 * zeta.^3 .* f_zeta + fp_zeta .* zeta4;
	decrs_dzeta = dfzeta4_dzeta .* gcrs - fp_zeta .* macrs * XC.fsec_inv;

	% Add LSD correlation functional to GGA exchange functional
	S.e_xc = S.e_xc + ecrs;
	vxcadd = ecrs - rs * XC.third .* decrs_drs - zeta .* decrs_dzeta;
	v_xc(:,1) = v_xc(:,1) + vxcadd + decrs_dzeta;
	v_xc(:,2) = v_xc(:,2) + vxcadd - decrs_dzeta;
	%[S.e_xc(1:20) v_xc(1:20,1) v_xc(1:20,2)] 
	% -----------------------------------------------------------------------------
	% Eventually add the GGA correlation part of the PBE functional

	% The definition of phi has been slightly changed, because
	% the original PBE one gives divergent behaviour for fully polarized points

	phi_zeta = ( zetpm1_3 .* (1.0 + zeta * XC.alpha_zeta) + zetmm1_3 .* (1.0 - zeta * XC.alpha_zeta)   ) * 0.5;
	phip_zeta = (zetpm1_3 - zetmm1_3) * XC.third * XC.alpha_zeta;
	phi_zeta_inv = 1.0./phi_zeta;
	phi_logder = phip_zeta .* phi_zeta_inv;
	phi3_zeta = phi_zeta .* phi_zeta .* phi_zeta;
	gamphi3inv = XC.gamma_inv * phi_zeta_inv .* phi_zeta_inv .* phi_zeta_inv;
	
	%        From ec to bb
	bb = ecrs .* gamphi3inv;
	dbb_drs = decrs_drs .* gamphi3inv;
	dbb_dzeta = gamphi3inv .* (decrs_dzeta - 3.0 * ecrs .* phi_logder);
	
	% From bb to cc
	exp_pbe = exp(-bb);
	cc = 1.0./(exp_pbe - 1.0);
	dcc_dbb = cc .* cc .* exp_pbe;
	dcc_drs = dcc_dbb .* dbb_drs;
	dcc_dzeta = dcc_dbb .* dbb_dzeta;

	% From cc to aa
	coeff_aa = XC.beta * XC.gamma_inv * phi_zeta_inv .* phi_zeta_inv;
	aa = coeff_aa .* cc;
	daa_drs = coeff_aa .* dcc_drs;
	daa_dzeta = -2.0 * aa .* phi_logder + coeff_aa .* dcc_dzeta;

	% Introduce tt : do not assume that the spin-dependent gradients are collinear
	grrho2 = sigma(:,1);
	dtt_dg = 2.0 * rhotot_inv .* rhotot_inv .* rhom1_3 * XC.coeff_tt;
	% Note that tt is (the t variable of PBE divided by phi) squared
	tt = 0.5 * grrho2 .* dtt_dg;

	% Get xx from aa and tt
	xx = aa .* tt;
	dxx_drs = daa_drs .* tt;
	dxx_dzeta = daa_dzeta .* tt;
	dxx_dtt = aa;

	% From xx to pade
	pade_den = 1.0./(1.0 + xx .* (1.0 + xx));
	pade = (1.0 + xx) .* pade_den;
	dpade_dxx = -xx .* (2.0 + xx) .* (pade_den.^2);
	dpade_drs = dpade_dxx .* dxx_drs;
	dpade_dtt = dpade_dxx .* dxx_dtt;
	dpade_dzeta = dpade_dxx .* dxx_dzeta;

	% From pade to qq
	coeff_qq = tt .* phi_zeta_inv .* phi_zeta_inv;
	qq = coeff_qq .* pade;
	dqq_drs = coeff_qq .* dpade_drs;
	dqq_dtt = pade .* phi_zeta_inv .* phi_zeta_inv + coeff_qq .* dpade_dtt;
	dqq_dzeta = coeff_qq .* (dpade_dzeta - 2.0 * pade .* phi_logder);

	% From qq to rr
	arg_rr = 1.0 + XC.beta * XC.gamma_inv * qq;
	div_rr = 1.0./arg_rr;
	rr = XC.gamma * log(arg_rr);
	drr_dqq = XC.beta * div_rr;
	drr_drs = drr_dqq .* dqq_drs;
	drr_dtt = drr_dqq .* dqq_dtt;
	drr_dzeta = drr_dqq .* dqq_dzeta;

	% From rr to hh
	hh = phi3_zeta .* rr;
	dhh_drs = phi3_zeta .* drr_drs;
	dhh_dtt = phi3_zeta .* drr_dtt;
	dhh_dzeta = phi3_zeta .* (drr_dzeta + 3.0 * rr .* phi_logder);
	% The GGA correlation energy is added
	S.e_xc = S.e_xc + hh;

	% From hh to the derivative of the energy wrt the density
	drhohh_drho = hh - XC.third * rs .* dhh_drs - zeta .* dhh_dzeta - (7.0/3.0) * tt .* dhh_dtt;
	v_xc(:,1) = v_xc(:,1) + drhohh_drho + dhh_dzeta;
	v_xc(:,2) = v_xc(:,2) + drhohh_drho - dhh_dzeta;
	
   
	% From hh to the derivative of the energy wrt to the gradient of the
	% density, divided by the gradient of the density
	% (The v3.3 definition includes the division by the norm of the gradient)

	S.dvxcdgrho(:,1) = rho(:,1) .* dtt_dg .* dhh_dtt;
	S.dvxcdgrho(:,2) = dvxcdgrho1(:,1) ;
	S.dvxcdgrho(:,3) = dvxcdgrho1(:,2) ;
	
	if S.cell_typ ~= 2
		Vxc_temp = S.grad_1 * (S.dvxcdgrho.*drho_1) + S.grad_2 * (S.dvxcdgrho.*drho_2) + S.grad_3 * (S.dvxcdgrho.*drho_3);
	else
		Vxc_temp = S.lapc_T(1,1)*S.grad_1*(S.dvxcdgrho.*drho_1) + S.lapc_T(2,2)*S.grad_2*(S.dvxcdgrho.*drho_2) + S.lapc_T(3,3)*S.grad_3*(S.dvxcdgrho.*drho_3) +...
				   S.lapc_T(2,1)*S.grad_1*(S.dvxcdgrho.*drho_2) + S.lapc_T(2,1)*S.grad_2*(S.dvxcdgrho.*drho_1) + S.lapc_T(3,2)*S.grad_2*(S.dvxcdgrho.*drho_3) +...
				   S.lapc_T(3,2)*S.grad_3*(S.dvxcdgrho.*drho_2) + S.lapc_T(3,1)*S.grad_1*(S.dvxcdgrho.*drho_3) + S.lapc_T(3,1)*S.grad_3*(S.dvxcdgrho.*drho_1) ;
	end
	S.Vxc(:,1) = v_xc(:,1) - Vxc_temp(:,2) - Vxc_temp(:,1);
	S.Vxc(:,2) = v_xc(:,2) - Vxc_temp(:,3) - Vxc_temp(:,1);             
end


