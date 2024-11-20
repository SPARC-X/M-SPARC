function [ec, vc, v2c, v3c] = rscanc(rho, sigma, tau)
% @file    scanx.m
% @brief   This file contains the functions computing the correlation energy density \epsilon_c and the potential V
% @authors Boqin Zhang <bzhang376@gatech.edu>
%          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% Reference:
% Bartók, Albert P., and Jonathan R. Yates. "Regularized SCAN functional." 
% The Journal of chemical physics 150, no. 16 (2019): 161101.
% Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
% ==============================================================================================
% rho is electron density n, nnr*1 vector 
% !!! attention: rho(rho < S.xc_rhotol) = S.xc_rhotol; process rho before
% calling the function!
% Drho is gradient of electron density n, nnr*3 vectors
% tau is kinetic energy density, nnr*1 vector
    normDrho = sigma.^0.5;
    [s, alphaP, DsDn, DsDDn, DalphaPDn, DalphaPDDn, DalphaPDtau] = basicrscanvariables(rho, normDrho, tau);
    [ec, vc, v2c, v3c] = correlationrSCAN(rho, s, alphaP, DsDn, DsDDn, DalphaPDn, DalphaPDDn, DalphaPDtau);
    v2c = v2c ./ normDrho; % D(rho*Ec)/D(|grad rho|) / |grad rho|
end
    
function [ec, vc, v2c, v3c] = correlationrSCAN(rho, s, alphaP, DsDn, DsDDn, DalphaPDn, DalphaPDDn, DalphaPDtau)
    zeta = 0; % now there is no spin!
    phi = ((1 + zeta).^(2/3) + (1 - zeta).^(2/3)) / 2; % now there is no spin!
    rs = (0.75./(pi*rho)).^(1/3) ;
    %% compute epsilon_c^0 when alpha \approx 0
    b1c = 0.0285764;
    b2c = 0.0889;
    b3c = 0.125541;
    ecLDA0 = -b1c ./ (1 + b2c*(rs.^0.5) + b3c*rs);
    dx = 0.5*((1 + zeta).^(4/3) + (1 - zeta).^(4/3)); % no spin, it should be 1
%     cx = -3/(4*pi)*(9*pi/4)^(1/3)*dx;
    cx0 = -3/(4*pi)*(9*pi/4)^(1/3);
    Gc = (1 - 2.3631*(dx - 1)) .* (1 - zeta^12);
    w0 = exp(-ecLDA0/b1c) - 1;
    betaConst = 0.06672455060314922;
%     betaRsInf = 0.066725*0.1/0.1778;
    betaRsInf = betaConst*0.1/0.1778;
    f0 = -0.9;
%     xiInf = (3*pi^2/16)^(2/3) * (betaRsInf*phi/(cx - f0)); % should be 0.128026 if zeta = 0
    xiInf0 = (3*pi^2/16)^(2/3) * (betaRsInf*1/(cx0 - f0));
%     gInf = (1 + 4*xiInf*s.^2).^(-0.25);
    gInf0s = (1 + 4*xiInf0*s.^2).^(-0.25);
    H0 = b1c*log(1 + w0.*(1 - gInf0s));
    ec0 = (ecLDA0 + H0).*Gc;
    %% compute epsilon_c^1 when alpha \approx 1
    sqr_rs = rs .^ 0.5;
    rsm1_2 = 1 ./ sqr_rs;
%     beta = 0.066725* (1 + 0.1*rs) ./ (1 + 0.1778*rs);
    beta = betaConst* (1 + 0.1*rs) ./ (1 + 0.1778*rs);
    % epsilon_c LSDA1
    % correlation parameters
	p = 1 ;
% 	AA = 0.031091 ;
	AA = 0.0310907 ;
	alpha1 = 0.21370 ;
	beta1 = 7.5957 ;
	beta2 = 3.5876 ;
	beta3 = 1.6382 ;
	beta4 = 0.49294 ;
    % conpute
    ec0_q0 = -2.0 * AA * (1.0 + alpha1*rs);
	ec0_q1 = 2.0 * AA *(beta1*sqr_rs + beta2*rs + beta3*rs.*sqr_rs + beta4*rs.*rs);
	ec0_q1p = AA * (beta1*rsm1_2 + 2.0*beta2 + 3.0*beta3*sqr_rs + 4.0*beta4*rs);
	ec0_den = 1.0 ./ (ec0_q1.*ec0_q1 + ec0_q1);
	ec0_log = -log(ec0_q1.*ec0_q1 .* ec0_den);
	ecrs0 = ec0_q0 .* ec0_log;
    ec_lsda1 = ecrs0;
	Dec_lsda1Drs = -2.0 * AA * alpha1 * ec0_log - ec0_q0 .* ec0_q1p .* ec0_den;
    
    r = 0.031091;
	w1 = exp(-ec_lsda1 ./ (r*phi.^3)) - 1;
    A = beta ./ (r*w1);
    t = (3*pi^2/16)^(1/3) * s./(phi*sqr_rs);
    g = (1 + 4*A.*t.*t).^(-0.25);
    H1 = r*phi.^3*log(1 + w1.*(1 - g));
    ec1 = ec_lsda1 + H1;
    %% interpolate and extrapolate epsilon_c
    c1c = 0.64;
    c2c = 1.5;
    dc = 0.7;
    alphaG25 = (alphaP > 2.5);
    alpha0To25 = ((alphaP >= 0.0) & (alphaP <= 2.5));
    alphaL0 = (alphaP < 0.0);
    fc = zeros(size(rho, 1), size(rho, 2));
    fc(alphaG25) = -dc*exp(c2c ./ (1 - alphaP(alphaG25)));
    fc(alpha0To25) = 1.0 + (-0.64)*alphaP(alpha0To25) + (-0.4352)*alphaP(alpha0To25).^2 + (-1.535685604549)*alphaP(alpha0To25).^3 ...
        + 3.061560252175*alphaP(alpha0To25).^4 + (-1.915710236206)*alphaP(alpha0To25).^5 + 0.516884468372*alphaP(alpha0To25).^6 ...
        + (-0.051848879792)*alphaP(alpha0To25).^7; % new formula for the interval
    fc(alphaL0) = exp(-c1c*alphaP(alphaL0) ./ (1 - alphaP(alphaL0)));
    ec = ec1 + fc.*(ec0 - ec1);
    %% compute variation of epsilon_c^0
    DzetaDn = 0; % no spin
    DrsDn = -4*pi/9*(4*pi/3*rho).^(-4/3);
    DdxDn = (4/3*(1 + zeta).^(1/3) - 4/3*(1 - zeta).^(1/3)).*DzetaDn; % when there is no spin, it should be 0
    DGcDn = -2.3631*DdxDn.*(1 - zeta.^12) + (1 - 2.3631*(dx - 1)).*(12*zeta.^11.*DzetaDn);
    DgInf0sDs = -0.25*(1 + 4*xiInf0*s.*s).^(-1.25) .* (4*xiInf0*2*s);
    DgInf0sDn = DgInf0sDs .* DsDn;
    DgInf0sDDn = DgInf0sDs .* DsDDn;
    DecLDA0Dn = b1c*(0.5*b2c*rs.^(-0.5) + b3c)./((1 + b2c*rs.^0.5 + b3c*rs).^2).*DrsDn;
    Dw0Dn = (w0 + 1) .* (-DecLDA0Dn/b1c);
    DH0Dn = b1c*(Dw0Dn.*(1 - gInf0s) - w0.*DgInf0sDn) ./ (1 + w0.*(1 - gInf0s));
    DH0DDn = b1c* (-w0.*DgInf0sDDn) ./ (1 + w0.*(1 - gInf0s));
    Dec0Dn = (DecLDA0Dn + DH0Dn).*Gc + (ecLDA0 + H0).*DGcDn;
    Dec0DDn = DH0DDn.*Gc;
    %% compute variation of epsilon_c^1
    Dec_lsda1Dn = - (rs./rho/3).*(-2*AA*alpha1*log(1+1./(2*AA*( beta1*sqr_rs + beta2*rs + beta3*(rs.^1.5) + beta4*(rs.^(p+1))))) ...
		- ((-2*AA*(1+alpha1*rs)).*(AA*( beta1*(rs.^-0.5)+ 2*beta2 + 3*beta3*(rs.^0.5) + 2*(p+1)*beta4*(rs.^p) ))) ...
		./((2*AA*( beta1*sqr_rs + beta2*rs + beta3*(rs.^1.5) + beta4*(rs.^(p+1)) ) ) ...
		.*(2*AA*( beta1*sqr_rs + beta2*rs + beta3*(rs.^1.5) + beta4*(rs.^(p+1)) ) )+(2*AA*( beta1*(rs.^0.5) + beta2*rs + beta3*(rs.^1.5) + beta4*(rs.^(p+1)) ) )) ) ;% from LDA_PW(S)
    DbetaDn = 0.066725*(0.1*(1+0.1778*rs) - 0.1778*(1+0.1*rs)) ./ ((1+0.1778*rs).^2) .* DrsDn;
    DphiDn = 0.5*(2/3*(1+zeta).^(-1/3) - 2/3*(1-zeta).^(-1/3)).*DzetaDn; % no spin, should be 0
    DtDn = (3*pi^2/16)^(1/3)*(phi.*sqr_rs.*DsDn - s.*(DphiDn.*sqr_rs + phi.*DrsDn./(2*sqr_rs))) ./ (phi.^2.*rs);
    DtDDn = t.*DsDDn./s;
    Dw1Dn = (w1 + 1) .* (-(r*phi.^3.*Dec_lsda1Dn - r*ec_lsda1.*(3.*phi.^2.*DphiDn)) ./ ((r*phi.^3).^2));
    DADn = (w1.*DbetaDn - beta.*Dw1Dn) ./ (r*w1.^2);
    DgDn = -0.25*(1 + 4*A.*t.*t).^(-1.25) .* (4*(DADn.*t.*t + 2*A.*t.*DtDn));
    DgDDn = -0.25*(1 + 4*A.*t.*t).^(-1.25) .* (4*2*A.*t.*DtDDn);
    DH1Dn = r*(3*phi.^2.*DphiDn.*log(1 + w1.*(1 - g)) + phi.^3.*(Dw1Dn.*(1 - g) - w1.*DgDn) ./ (1 + w1.*(1 - g)));
    DH1DDn = r* (phi.^3.*(-w1.*DgDDn) ./ (1 + w1.*(1 - g)));
    Dec1Dn = Dec_lsda1Dn + DH1Dn;
    Dec1DDn = DH1DDn;
    %% variant of fc and ec
    DfcDalpha = zeros(size(rho, 1), size(rho, 2));
    DfcDalpha(alphaG25) = fc(alphaG25).*(c2c./(1 - alphaP(alphaG25)).^2);
    DfcDalpha(alpha0To25) = (-0.64) + (-0.4352)*alphaP(alpha0To25)*2 + (-1.535685604549)*alphaP(alpha0To25).^2*3 ...
        + 3.061560252175*alphaP(alpha0To25).^3*4 + (-1.915710236206)*alphaP(alpha0To25).^4*5 + 0.516884468372*alphaP(alpha0To25).^5*6 ...
        + (-0.051848879792)*alphaP(alpha0To25).^6*7; % new formula for the interval;
    DfcDalpha(alphaL0) = fc(alphaL0).*(-c1c./(1 - alphaP(alphaL0)).^2);
    DfcDn = DfcDalpha.*DalphaPDn;
    DfcDDn = DfcDalpha.*DalphaPDDn;
    DfcDtau = DfcDalpha.*DalphaPDtau;
    DepsiloncDn = Dec1Dn + fc.*(Dec0Dn -Dec1Dn) + DfcDn.*(ec0 - ec1);
    DepsiloncDDn = Dec1DDn + fc.*(Dec0DDn -Dec1DDn) + DfcDDn.*(ec0 - ec1);
    DepsiloncDtau = DfcDtau.*(ec0 - ec1);
    vc = ec + rho.*DepsiloncDn;
    v2c = rho.*DepsiloncDDn; % D(rho*Ec)/D(|grad rho|)
    v3c = rho.*DepsiloncDtau;
end

function [s, alphaP, DsDn, DsDDn, DalphaPDn, DalphaPDDn, DalphaPDtau] = basicrscanvariables(rho, normDrho, tau)
    %% variables
    s = normDrho ./ (2*(3*pi^2)^(1/3)*rho.^(4/3));
    tauw = normDrho.^2 ./ (8*rho);
    tauUnif = 3/10*(3*pi^2)^(2/3) * rho.^(5/3);
    alpha = (tau - tauw) ./ (tauUnif + 1e-4); % when there is no spin, ds = 1
    alphaP = alpha.^3 ./ (alpha.^2 + 1e-3);
    %% variation
    DsDn = -2*normDrho ./ (3*(3*pi^2)^(1/3).*rho.^(7/3));
    DsDDn = 1 ./ (2*(3*pi^2)^(1/3)*rho.^(4/3));
    DtauwDn = -normDrho.^2 ./ (8*rho.^2);
    DtauwDDn = normDrho ./ (4*rho);
    DtauUnifDn = (3*pi^2)^(2/3)/2 * rho.^(2/3);
    DalphaDn = (-DtauwDn.*(tauUnif + 1e-4) - (tau - tauw).*DtauUnifDn) ./ (tauUnif + 1e-4).^2;
    DalphaDDn = (-DtauwDDn) ./ (tauUnif + 1e-4);
    DalphaDtau = 1 ./ (tauUnif + 1e-4);
    DalphaPDalpha = (3*alpha.^2.*(alpha.^2 + 1e-3) - alpha.^3.*(2*alpha)) ./ (alpha.^2 + 1e-3).^2;
    DalphaPDn = DalphaPDalpha.*DalphaDn;
    DalphaPDDn = DalphaPDalpha.*DalphaDDn;
    DalphaPDtau = DalphaPDalpha.*DalphaDtau;
end
