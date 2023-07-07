function [ec,vc,v2c,v3c] = r2scanc_spin(rho,sigma,tau)
% @file    r2scanc_spin.m
% @brief   This file contains the functions computing the spin-polarized correlation energy density \epsilon_c and the potential V
% @authors Boqin Zhang <bzhang376@gatech.edu>
%          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% Reference:
% Furness, James W., Aaron D. Kaplan, Jinliang Ning, John P. Perdew, and Jianwei Sun. 
% "Accurate and numerically efficient r2SCAN meta-generalized gradient approximation." 
% The journal of physical chemistry letters 11, no. 19 (2020): 8208-8215.
% Physical review letters 115, no. 3 (2015): 036402.
% Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
% ==============================================================================================
    normDrho = sigma.^0.5; % col1: grad of n; col2: grad of n up; col3: grad of n dn
    [s, p, alpha, zeta, DsDn, DsDDn, DpDn, DpDDn, DalphaDnup, DalphaDndn, DalphaDDn, DalphaDtau, DzetaDnup, DzetaDndn] ...
        = co_basicr2SCANvariables(rho, normDrho(:, 1), tau(:, 1));
    [ec, vc, v2c, v3c] = correlationSr2SCAN(rho(:, 1), s, p, alpha, zeta, DsDn, DsDDn, DpDn, DpDDn, DalphaDnup, DalphaDndn, DalphaDDn, DalphaDtau, DzetaDnup, DzetaDndn);
    v2c = v2c./normDrho(:, 1);
    v3c = [v3c, v3c];
end

function [s, p, alpha, zeta, DsDn, DsDDn, DpDn, DpDDn, DalphaDnup, DalphaDndn, DalphaDDn, DalphaDtau, DzetaDnup, DzetaDndn] ...
    = co_basicr2SCANvariables(rho, normDrho, tau)
%     normDrho = (Drho(:,1).^2 + Drho(:,2).^2 + Drho(:,3).^2).^0.5;
    % rho at here has 3 cols
    %% value
    s = normDrho ./ (2*(3*pi^2)^(1/3)*rho(:, 1).^(4/3));
    p = s .^2;
    zeta = (rho(:, 2) - rho(:, 3)) ./ rho(:, 1);
    DzetaDnup = 2*rho(:, 3) ./ (rho(:, 1).^2);
    DzetaDndn = -2*rho(:, 2) ./ (rho(:, 1).^2);
    tauw = normDrho.^2 ./ (8*rho(:, 1));
    ds = ((1+zeta).^(5/3) + (1-zeta).^(5/3)) / 2;
    tauUnif = 3/10*(3*pi^2)^(2/3) * rho(:, 1).^(5/3) .* ds;
    eta = 0.001;
    alpha = (tau - tauw) ./ (tauUnif + eta*tauw); % when there is no spin, ds = 1
    %% variation
    DsDn = -2*normDrho ./ (3*(3*pi^2)^(1/3).*rho(:, 1).^(7/3));
    DpDn = 2*s .* DsDn;
    DsDDn = 1 ./ (2*(3*pi^2)^(1/3)*rho(:, 1).^(4/3));
    DpDDn = 2*s .* DsDDn;
    DtauwDn = -normDrho.^2 ./ (8*rho(:, 1).^2);
    DtauwDDn = normDrho ./ (4*rho(:, 1));
%     DtauUnifDn = (3*pi^2)^(2/3)/2 * rho.^(2/3);
    DdsDnup = 5/3 * ((1+zeta).^(2/3) - (1-zeta).^(2/3)) .* DzetaDnup / 2;
    DdsDndn = 5/3 * ((1+zeta).^(2/3) - (1-zeta).^(2/3)) .* DzetaDndn / 2;
    DtauUnifDnup = (3*pi^2)^(2/3)/2 * rho(:, 1).^(2/3) .* ds + 3/10*(3*pi^2)^(2/3) * rho(:, 1).^(5/3) .* DdsDnup;
    DtauUnifDndn = (3*pi^2)^(2/3)/2 * rho(:, 1).^(2/3) .* ds + 3/10*(3*pi^2)^(2/3) * rho(:, 1).^(5/3) .* DdsDndn;
%     DalphaDn = (-DtauwDn.*(tauUnif + eta*tauw) - (tau - tauw).*(DtauUnifDn + eta*DtauwDn)) ./ ((tauUnif + eta*tauw).^2);
    DalphaDnup = (-DtauwDn.*(tauUnif + eta*tauw) - (tau - tauw).*(DtauUnifDnup + eta*DtauwDn)) ./ ((tauUnif + eta*tauw).^2);
    DalphaDndn = (-DtauwDn.*(tauUnif + eta*tauw) - (tau - tauw).*(DtauUnifDndn + eta*DtauwDn)) ./ ((tauUnif + eta*tauw).^2);
    DalphaDDn = (-DtauwDDn.*(tauUnif + eta*tauw) - (tau - tauw).*eta.*DtauwDDn) ./ ((tauUnif + eta*tauw).^2);
    DalphaDtau = 1 ./ (tauUnif + eta*tauw);
end

function [ec, Vc1, Vc2, Vc3] = correlationSr2SCAN(rho, s, p, alpha, zeta, DsDn, DsDDn, DpDn, DpDDn, DalphaDnup, DalphaDndn, DalphaDDn, DalphaDtau, DzetaDnup, DzetaDndn)
% all input variables of this function only have one column. 
% ec has only one column; Vc1 has two columns; Vc2 and Vc3 have only one column.
    phi = ((1 + zeta).^(2/3) + (1 - zeta).^(2/3)) / 2;
    
    rs = (0.75./(pi*rho)).^(1/3) ;
    %% compute epsilon_c^0 when alpha \approx 0
    b1c = 0.0285764;
    b2c = 0.0889;
    b3c = 0.125541;
    ecLDA0 = -b1c ./ (1 + b2c*(rs.^0.5) + b3c*rs);
    dx = ((1 + zeta).^(4/3) + (1 - zeta).^(4/3)) / 2; % no spin, it should be 1
    cx0 = -3/(4*pi)*(9*pi/4)^(1/3); % unused variable
    Gc = (1 - 2.3631*(dx - 1)) .* (1 - zeta.^12);
    w0 = exp(-ecLDA0/b1c) - 1;
    betaConst = 0.06672455060314922;
    betaRsInf = betaConst*0.1/0.1778; % the value is estimated
    f0 = -0.9;
    chiInf = (3*pi^2/16)^(2/3) * (betaRsInf*1/(cx0 - f0)); % the value is estimated, % should be 0.128026 if zeta = 0
%     chiInf = 0.128025; % there is a more accurate formula (seems like the previous formula). It is an estimate value.
%     gInf0s = (1 + 4*xiInf0.*(s.^2)).^(-0.25); % different formula, to be replaced.
    gInf0s = (1 + 4*chiInf.*(s.^2)).^(-0.25);
    H0 = b1c*log(1 + w0.*(1 - gInf0s));
    ec0 = (ecLDA0 + H0).*Gc;
    %% compute epsilon_c^1 when alpha \approx 1
    sqr_rs = rs .^ 0.5;
    rsm1_2 = 1 ./ sqr_rs;
%     beta = 0.066725* (1 + 0.1*rs) ./ (1 + 0.1778*rs);
    beta = betaConst* (1 + 0.1*rs) ./ (1 + 0.1778*rs);
    % epsilon_c LSDA1
    % correlation parameters
% 	p = 1 ;
% 	AA = 0.031091 ;
    AAec0 = 0.0310907 ;
    alpha1ec0 = 0.21370 ;
    beta1ec0 = 7.5957 ;
    beta2ec0 = 3.5876 ;
    beta3ec0 = 1.6382 ;
    beta4ec0 = 0.49294 ;
    % conpute
    ec0_q0 = -2.0 * AAec0 * (1.0 + alpha1ec0*rs);
    ec0_q1 = 2.0 * AAec0 *(beta1ec0*sqr_rs + beta2ec0*rs + beta3ec0*rs.*sqr_rs + beta4ec0*rs.*rs);
    ec0_q1p = AAec0 * (beta1ec0*rsm1_2 + 2.0*beta2ec0 + 3.0*beta3ec0*sqr_rs + 4.0*beta4ec0*rs);
    ec0_den = 1.0 ./ (ec0_q1.*ec0_q1 + ec0_q1);
    ec0_log = -log(ec0_q1.*ec0_q1 .* ec0_den);
    ecrs0 = ec0_q0 .* ec0_log;
    decrs0_drs = -2.0 * AAec0 * alpha1ec0 * ec0_log - ec0_q0 .* ec0_q1p .* ec0_den;
    
    AAmac = 0.0168869 ;
    alpha1mac = 0.11125 ;
    beta1mac = 10.357 ;
    beta2mac = 3.6231 ;
    beta3mac = 0.88026 ;
    beta4mac = 0.49671 ;
    mac_q0 = -2.0 * AAmac * (1.0 + alpha1mac * rs);
    mac_q1 = 2.0 * AAmac * (beta1mac * sqr_rs + beta2mac * rs + beta3mac * rs .* sqr_rs + beta4mac * rs .* rs);
    mac_q1p = AAmac * (beta1mac * rsm1_2 + 2 * beta2mac + 3 * beta3mac * sqr_rs + 4 * beta4mac * rs);
    mac_den = 1.0./(mac_q1 .* mac_q1 + mac_q1);
    mac_log = -log( mac_q1 .* mac_q1 .* mac_den );
    macrs = mac_q0 .* mac_log;
    dmacrs_drs = -2.0 * AAmac * alpha1mac * mac_log - mac_q0 .* mac_q1p .* mac_den;
    
    AAec1 = 0.01554535 ;
    alpha1ec1 = 0.20548 ;
    beta1ec1 = 14.1189 ;
    beta2ec1 = 6.1977 ;
    beta3ec1 = 3.3662 ;
    beta4ec1 = 0.62517 ;
    ec1_q0 = -2.0 * AAec1 * (1.0 + alpha1ec1 * rs);
    ec1_q1 = 2.0 * AAec1 * (beta1ec1 * sqr_rs + beta2ec1 * rs + beta3ec1 * rs .* sqr_rs + beta4ec1 * rs .* rs);
    ec1_q1p = AAec1 * (beta1ec1 * rsm1_2 + 2 * beta2ec1 + 3 * beta3ec1 * sqr_rs + 4 * beta4ec1 * rs);
    ec1_den = 1.0./(ec1_q1 .* ec1_q1 + ec1_q1);
    ec1_log = -log( ec1_q1 .* ec1_q1 .* ec1_den );
    ecrs1 = ec1_q0 .* ec1_log;
    decrs1_drs = -2.0 * AAec1 * alpha1ec1 * ec1_log - ec1_q0 .* ec1_q1p .* ec1_den;
    
    
    f_zeta = ((1.0+zeta).^(4/3) + (1.0-zeta).^(4/3) - 2.0 ) / (2^(4/3) - 2);
    fp_zeta = ((1.0+zeta).^(1/3) - (1.0-zeta).^(1/3)) * 4.0/3.0 / (2^(4/3) - 2);
    zeta4 = zeta.^4;

    gcrs = ecrs1 - ecrs0 + macrs/1.709921;
    ec_lsda1 = ecrs0 + f_zeta .* (zeta4 .* gcrs - macrs/1.709921);
    declsda1_drs = decrs0_drs + f_zeta .* (zeta4 .* (decrs1_drs - decrs0_drs + dmacrs_drs/1.709921) - dmacrs_drs/1.709921);
    
    r = 0.031091;
    w1 = exp(-ec_lsda1 ./ (r*phi.^3)) - 1;
%     A = beta ./ (r*w1); % the variable is replaced by y below
    t = (3*pi^2/16)^(1/3) * s./(phi.*sqr_rs);
%     g = (1 + 4*A.*t.*t).^(-0.25); % the formula is replaced by a new formula

    y = beta./(r*w1) .* (t.^2); % new variable, y = A*t*t
    deltafc2 = 1*(-0.64) + 2*(-0.4352) + 3*(-1.535685604549) + 4*3.061560252175 + 5*(-1.915710236206) + 6*0.516884468372 + 7*(-0.051848879792); 
    % new variable, a constant (or a polynomial?)
    ds = ((1 + zeta).^(5/3) + (1 - zeta).^(5/3)) / 2; % new variable. no spin, it should be 1
    eta = 0.001;
    dp2 = 0.361;
    ec_lsda0 = ecLDA0 .* Gc; % new variable, decLDA0drs is in the formula of DecLDA0Dn
    declsda0_drs = b1c*(0.5*b2c*rs.^(-0.5) + b3c)./((1 + b2c*rs.^0.5 + b3c*rs).^2) .* Gc; % new variable
    
    deltayPart1 = deltafc2 ./ (27*r*ds.*(phi.^3).*w1);
    deltayPart2 = 20*rs.*(declsda0_drs - declsda1_drs) - 45*eta*(ec_lsda0 - ec_lsda1);
    deltayPart3 = p.*exp(-p.^2/dp2^4); % new variable
    
    deltay = deltayPart1 .* deltayPart2 .* deltayPart3;
    
    g = (1 + 4*(y - deltay)).^(-0.25); % new formula
    
    H1 = r*phi.^3.*log(1 + w1.*(1 - g));
    ec1 = ec_lsda1 + H1; % ec_lsda1 in scan is ec_lsda in r2scan
    %% interpolate and extrapolate epsilon_c
    c1c = 0.64;
    c2c = 1.5;
    dc = 0.7;
%     alphaG1 = (alpha > 1);
%     alphaE1 = (alpha == 1);
%     alphaL1 = (alpha < 1); replaced by new formula
    alphaG25 = (alpha > 2.5);
    alpha0To25 = ((alpha >= 0.0) & (alpha <= 2.5));
    alphaL0 = (alpha < 0.0);
    fc = zeros(size(rho, 1), size(rho, 2));
    fc(alphaG25) = -dc*exp(c2c ./ (1 - alpha(alphaG25)));
    fc(alpha0To25) = 1.0 + (-0.64)*alpha(alpha0To25) + (-0.4352)*alpha(alpha0To25).^2 + (-1.535685604549)*alpha(alpha0To25).^3 ...
        + 3.061560252175*alpha(alpha0To25).^4 + (-1.915710236206)*alpha(alpha0To25).^5 + 0.516884468372*alpha(alpha0To25).^6 ...
        + (-0.051848879792)*alpha(alpha0To25).^7; % new formula for the interval
    fc(alphaL0) = exp(-c1c*alpha(alphaL0) ./ (1 - alpha(alphaL0)));
    ec = ec1 + fc.*(ec0 - ec1);
    %% compute variation of epsilon_c^0
%     DzetaDn = 0; % no spin
    DrsDn = -4*pi/9*(4*pi/3*rho).^(-4/3);
    DdxDnup = (4/3*(1 + zeta).^(1/3) - 4/3*(1 - zeta).^(1/3))/2.*DzetaDnup; 
    DdxDndn = (4/3*(1 + zeta).^(1/3) - 4/3*(1 - zeta).^(1/3))/2.*DzetaDndn; 
    DGcDnup = -2.3631*DdxDnup.*(1 - zeta.^12) + (1 - 2.3631*(dx - 1)).*(-12*zeta.^11.*DzetaDnup);
    DGcDndn = -2.3631*DdxDndn.*(1 - zeta.^12) + (1 - 2.3631*(dx - 1)).*(-12*zeta.^11.*DzetaDndn);
    DgInf0sDs = -0.25*(1 + 4*chiInf*s.*s).^(-1.25) .* (4*chiInf*2*s);
    DgInf0sDn = DgInf0sDs .* DsDn;
    DgInf0sDDn = DgInf0sDs .* DsDDn;
    DecLDA0Dn = b1c*(0.5*b2c*rs.^(-0.5) + b3c)./((1 + b2c*rs.^0.5 + b3c*rs).^2).*DrsDn;
    Dw0Dn = (w0 + 1) .* (-DecLDA0Dn/b1c);
    DH0Dn = b1c*(Dw0Dn.*(1 - gInf0s) - w0.*DgInf0sDn) ./ (1 + w0.*(1 - gInf0s));
    DH0DDn = b1c* (-w0.*DgInf0sDDn) ./ (1 + w0.*(1 - gInf0s));
%     Dec0Dn = (DecLDA0Dn + DH0Dn).*Gc + (ec0 + H0).*DGcDn;
    Dec0Dnup = (DecLDA0Dn + DH0Dn).*Gc + (ecLDA0 + H0).*DGcDnup;
    Dec0Dndn = (DecLDA0Dn + DH0Dn).*Gc + (ecLDA0 + H0).*DGcDndn;
    Dec0DDn = DH0DDn.*Gc;
    %% compute variation of epsilon_c^1
    dgcrs_drs = decrs1_drs - decrs0_drs + dmacrs_drs/1.709921;
    dec_lsda1_drs = decrs0_drs + f_zeta .* (zeta4 .* dgcrs_drs - dmacrs_drs/1.709921);
    DdsDnup = (5/3*(1 + zeta).^(2/3) - 5/3*(1 - zeta).^(2/3))/2.*DzetaDnup; 
    DdsDndn = (5/3*(1 + zeta).^(2/3) - 5/3*(1 - zeta).^(2/3))/2.*DzetaDndn; 
    
    dec0log_drs = -ec0_q1p .* ec0_den;
    dec0q0_drs = -2.0 * AAec0 * alpha1ec0;
    dec0q1p_drs = AAec0 * ((-0.5)*beta1ec0*rsm1_2./rs + 3.0*0.5*beta3ec0*rsm1_2 + 4.0*beta4ec0);
    dec0den_drs = -(2*ec0_q1.*ec0_q1p + ec0_q1p) ./ (ec0_q1.*ec0_q1 + ec0_q1).^2;
    d2ecrs0_drs2 = -2.0 * AAec0 * alpha1ec0 * dec0log_drs - (dec0q0_drs .* ec0_q1p .* ec0_den + ...
        ec0_q0 .* dec0q1p_drs .* ec0_den + ec0_q0 .* ec0_q1p .* dec0den_drs);
    
    dmaclog_drs = -mac_q1p .* mac_den;
    dmacq0_drs = -2.0 * AAmac * alpha1mac;
    dmacq1p_drs = AAmac * ((-0.5)*beta1mac*rsm1_2./rs + 3.0*0.5*beta3mac*rsm1_2 + 4.0*beta4mac);
    dmacden_drs = -(2*mac_q1.*mac_q1p + mac_q1p) ./ (mac_q1.*mac_q1 + mac_q1).^2;
    d2macrs_drs2 = -2.0 * AAmac * alpha1mac * dmaclog_drs - (dmacq0_drs .* mac_q1p .* mac_den + ...
        mac_q0 .* dmacq1p_drs .* mac_den + mac_q0 .* mac_q1p .* dmacden_drs);
    
    dec1log_drs = -ec1_q1p .* ec1_den;
    dec1q0_drs = -2.0 * AAec1 * alpha1ec1;
    dec1q1p_drs = AAec1 * ((-0.5)*beta1ec1*rsm1_2./rs + 3.0*0.5*beta3ec1*rsm1_2 + 4.0*beta4ec1);
    dec1den_drs = -(2*ec1_q1.*ec1_q1p + ec1_q1p) ./ (ec1_q1.*ec1_q1 + ec1_q1).^2;
    d2ecrs1_drs2 = -2.0 * AAec1 * alpha1ec1 * dec1log_drs - (dec1q0_drs .* ec1_q1p .* ec1_den + ...
        ec1_q0 .* dec1q1p_drs .* ec1_den + ec1_q0 .* ec1_q1p .* dec1den_drs);
    
    d2eclsda1_drs2 = d2ecrs0_drs2 + f_zeta .* (zeta4 .* (d2ecrs1_drs2 - d2ecrs0_drs2 + d2macrs_drs2/1.709921) - d2macrs_drs2/1.709921);
    d2eclsda1_drsdzeta = fp_zeta.* (zeta4 .* (decrs1_drs - decrs0_drs + dmacrs_drs/1.709921) - dmacrs_drs/1.709921) ...
        + f_zeta.*(4*zeta.^3 .* (decrs1_drs - decrs0_drs + dmacrs_drs/1.709921));
    Ddeclsda1_drsDnup = d2eclsda1_drs2.*DrsDn + d2eclsda1_drsdzeta.*DzetaDnup;
    Ddeclsda1_drsDndn = d2eclsda1_drs2.*DrsDn + d2eclsda1_drsdzeta.*DzetaDndn;
    % new variable, derives from decrs0_drs, dmacrs_drs, decrs1_drs formula above
    
    dfzeta4_dzeta = 4.0 * zeta.^3 .* f_zeta + fp_zeta .* zeta4;
    %     declsda1_dzeta = fp_zeta .* (zeta4 .* gcrs - macrs/1.709921) + f_zeta .* (4*zeta.^3.*gcrs);
    dec_lsda1_dzeta = dfzeta4_dzeta .* gcrs - fp_zeta .* macrs/1.709921;
%     vxcadd = ec_lsda1 - rs * 1/3 .* dec_lsda1_drs - zeta .* dec_lsda1_dzeta;
    Dec_lsda1Dnup = (- rs * 1/3 .* dec_lsda1_drs - zeta .* dec_lsda1_dzeta + dec_lsda1_dzeta) ./ rho(:, 1);
    Dec_lsda1Dndn = (- rs * 1/3 .* dec_lsda1_drs - zeta .* dec_lsda1_dzeta - dec_lsda1_dzeta) ./ rho(:, 1); % from LDA_PW(S)
%     fprintf("ecLSDA %10.8f, Vc1LSDA %10.8f %10.8f\n", ec_lsda1, Dec_lsda1Dnup, Dec_lsda1Dndn);
    DbetaDn = 0.066725*(0.1*(1+0.1778*rs) - 0.1778*(1+0.1*rs)) ./ ((1+0.1778*rs).^2) .* DrsDn;
%     DphiDn = 0.5*(2/3*(1+zeta).^(-1/3) - 2/3*(1-zeta).^(-1/3)).*DzetaDn; % no spin, should be 0
    DphiDnup = 0.5*(2/3*(1+zeta).^(-1/3) - 2/3*(1-zeta).^(-1/3)).*DzetaDnup;
    DphiDndn = 0.5*(2/3*(1+zeta).^(-1/3) - 2/3*(1-zeta).^(-1/3)).*DzetaDndn;
%     DtDn = (3*pi^2/16)^(1/3)*(phi.*sqr_rs.*DsDn - s.*(DphiDn.*sqr_rs + phi.*DrsDn./(2*sqr_rs))) ./ (phi.^2.*rs);
    DtDnup = (3*pi^2/16)^(1/3)*(phi.*sqr_rs.*DsDn - s.*(DphiDnup.*sqr_rs + phi.*DrsDn./(2*sqr_rs))) ./ (phi.^2.*rs);
    DtDndn = (3*pi^2/16)^(1/3)*(phi.*sqr_rs.*DsDn - s.*(DphiDndn.*sqr_rs + phi.*DrsDn./(2*sqr_rs))) ./ (phi.^2.*rs);
    DtDDn = t.*DsDDn./s;
%     Dw1Dn = (w1 + 1) .* (-(r*phi.^3.*Dec_lsda1Dn - r*ec_lsda1.*(3.*phi.^2.*DphiDn)) ./ ((r*phi.^3).^2));
%     dw1drs = (w1 + 1).*(-declsda1_drs ./ (r*phi.^3));
%     Dw1Dzeta = (w1 + 1).*(-declsda1_dzeta ./ (r*phi.^3) - ec_lsda1.*(-3)./ (r*phi.^4));
    Dw1Dnup = (w1 + 1) .* (-(r*phi.^3.*Dec_lsda1Dnup - r*ec_lsda1.*(3.*phi.^2.*DphiDnup)) ./ ((r*phi.^3).^2));
    Dw1Dndn = (w1 + 1) .* (-(r*phi.^3.*Dec_lsda1Dndn - r*ec_lsda1.*(3.*phi.^2.*DphiDndn)) ./ ((r*phi.^3).^2));
% %     DADn = (w1.*DbetaDn - beta.*Dw1Dn) ./ (r*w1.^2);
%     DADnup = (w1.*DbetaDn - beta.*Dw1Dnup) ./ (r*w1.^2);
%     DADndn = (w1.*DbetaDn - beta.*Dw1Dndn) ./ (r*w1.^2);
    DyDnup = (w1.*DbetaDn - beta.*Dw1Dnup) ./ (r*w1.^2) .* (t.^2) + beta./(r*w1) .* (2*t).*DtDnup; % new variable
    DyDndn = (w1.*DbetaDn - beta.*Dw1Dndn) ./ (r*w1.^2) .* (t.^2) + beta./(r*w1) .* (2*t).*DtDndn;
    DyDDn = beta./(r*w1) .* (2*t).*DtDDn; % new variable
    
    Declsda0Dnup = declsda0_drs.*DrsDn + ecLDA0.*DGcDnup;
    Declsda0Dndn = declsda0_drs.*DrsDn + ecLDA0.*DGcDndn;
    d2eclsda0_drs2 = b1c*((0.5*b2c*(-0.5)*rs.^(-1.5)).*((1 + b2c*rs.^0.5 + b3c*rs).^2) - ...
        (0.5*b2c*rs.^(-0.5) + b3c).*2.*(1 + b2c*rs.^0.5 + b3c*rs).*(0.5*b2c*rs.^(-0.5) + b3c)) ...
        ./ (((1 + b2c*rs.^0.5 + b3c*rs).^2).^2) .* Gc;
    d2eclsda0_drsdGc = b1c*(0.5*b2c*rs.^(-0.5) + b3c)./((1 + b2c*rs.^0.5 + b3c*rs).^2);
    Ddeclsda0_drsDnup = d2eclsda0_drs2.*DrsDn + d2eclsda0_drsdGc.*DGcDnup;
    Ddeclsda0_drsDndn = d2eclsda0_drs2.*DrsDn + d2eclsda0_drsdGc.*DGcDndn;

    d_deltayPart1_dnup =  deltafc2 ./ (27*r*(phi.^3).*w1) .* (-1) .*ds.^(-2) .* DdsDnup ...
        + deltafc2 ./ (27*r*ds.*w1) .* (-3).*phi.^(-4) .* DphiDnup ...
        + deltafc2 ./ (27*r*ds.*(phi.^3)) .* (-1).*w1.^(-2) .* Dw1Dnup;
    d_deltayPart1_dndn =  deltafc2 ./ (27*r*(phi.^3).*w1) .* (-1) .*ds.^(-2) .* DdsDndn ...
        + deltafc2 ./ (27*r*ds.*w1) .* (-3).*phi.^(-4) .* DphiDndn ...
        + deltafc2 ./ (27*r*ds.*(phi.^3)) .* (-1).*w1.^(-2) .* Dw1Dndn;

    d_deltayPart2_dnup = 20*(declsda0_drs - declsda1_drs).*DrsDn ...
        + 20*rs.*(Ddeclsda0_drsDnup - Ddeclsda1_drsDnup) ...
        - 45*eta*(Declsda0Dnup - Dec_lsda1Dnup);
    d_deltayPart2_dndn = 20*(declsda0_drs - declsda1_drs).*DrsDn ...
        + 20*rs.*(Ddeclsda0_drsDndn - Ddeclsda1_drsDndn) ...
        - 45*eta*(Declsda0Dndn - Dec_lsda1Dndn);
    d_deltayPart3_dp = exp(-p.^2/dp2^4) + p.*exp(-p.^2/dp2^4).*(-2*p/dp2^4);
    
    DdeltayDnup = d_deltayPart1_dnup.*deltayPart2.*deltayPart3 ...
        + deltayPart1.*d_deltayPart2_dnup.*deltayPart3 ...
        + deltayPart1.*deltayPart2.*d_deltayPart3_dp.*DpDn; % new variable
    DdeltayDndn = d_deltayPart1_dndn.*deltayPart2.*deltayPart3 ...
        + deltayPart1.*d_deltayPart2_dndn.*deltayPart3 ...
        + deltayPart1.*deltayPart2.*d_deltayPart3_dp.*DpDn; % new variable
    DdeltayDDn = deltayPart1.*deltayPart2.*d_deltayPart3_dp.*DpDDn;
%     DgDn = -0.25*(1 + 4*A.*t.*t).^(-1.25) .* (4*(DADn.*t.*t + 2*A.*t.*DtDn));
%     DgDnup = -0.25*(1 + 4*A.*t.*t).^(-1.25) .* (4*(DADnup.*t.*t + 2*A.*t.*DtDnup));
%     DgDndn = -0.25*(1 + 4*A.*t.*t).^(-1.25) .* (4*(DADndn.*t.*t + 2*A.*t.*DtDndn));
%     DgDDn = -0.25*(1 + 4*A.*t.*t).^(-1.25) .* (4*2*A.*t.*DtDDn);
    DgDnup = -0.25*(1 + 4*(y - deltay)).^(-1.25) .* (4*(DyDnup - DdeltayDnup)); % new formula
    DgDndn = -0.25*(1 + 4*(y - deltay)).^(-1.25) .* (4*(DyDndn - DdeltayDndn));
    DgDDn = -0.25*(1 + 4*(y - deltay)).^(-1.25) .* (4*(DyDDn - DdeltayDDn));
%     DH1Dn = r*(3*phi.^2.*DphiDn.*log(1 + w1.*(1 - g)) + phi.^3.*(Dw1Dn.*(1 - g) - w1.*DgDn) ./ (1 + w1.*(1 - g)));
    DH1Dnup = r*(3*phi.^2.*DphiDnup.*log(1 + w1.*(1 - g)) + phi.^3.*(Dw1Dnup.*(1 - g) - w1.*DgDnup) ./ (1 + w1.*(1 - g)));
    DH1Dndn = r*(3*phi.^2.*DphiDndn.*log(1 + w1.*(1 - g)) + phi.^3.*(Dw1Dndn.*(1 - g) - w1.*DgDndn) ./ (1 + w1.*(1 - g)));
    DH1DDn = r* (phi.^3.*(-w1.*DgDDn) ./ (1 + w1.*(1 - g)));
%     Dec1Dn = Dec_lsda1Dn + DH1Dn;
    Dec1Dnup = Dec_lsda1Dnup + DH1Dnup;
    Dec1Dndn = Dec_lsda1Dndn + DH1Dndn;
    Dec1DDn = DH1DDn;
    %% variant of fc and ec
    DfcDalpha = zeros(size(rho, 1), size(rho, 2));
    DfcDalpha(alphaG25) = fc(alphaG25).*(c2c./(1 - alpha(alphaG25)).^2);
    DfcDalpha(alpha0To25) = (-0.64) + (-0.4352)*alpha(alpha0To25)*2 + (-1.535685604549)*alpha(alpha0To25).^2*3 ...
        + 3.061560252175*alpha(alpha0To25).^3*4 + (-1.915710236206)*alpha(alpha0To25).^4*5 + 0.516884468372*alpha(alpha0To25).^5*6 ...
        + (-0.051848879792)*alpha(alpha0To25).^6*7; % new formula for the interval;
    DfcDalpha(alphaL0) = fc(alphaL0).*(-c1c./(1 - alpha(alphaL0)).^2);
    DfcDnup = DfcDalpha.*DalphaDnup;
    DfcDndn = DfcDalpha.*DalphaDndn;
    DfcDDn = DfcDalpha.*DalphaDDn;
    DfcDtau = DfcDalpha.*DalphaDtau;
%     DepsiloncDn = Dec1Dn + fc.*(Dec0Dn -Dec1Dn) + DfcDn.*(ec0 - ec1);
    DepsiloncDnup = Dec1Dnup + fc.*(Dec0Dnup -Dec1Dnup) + DfcDnup.*(ec0 - ec1);
    DepsiloncDndn = Dec1Dndn + fc.*(Dec0Dndn -Dec1Dndn) + DfcDndn.*(ec0 - ec1);
    DepsiloncDDn = Dec1DDn + fc.*(Dec0DDn -Dec1DDn) + DfcDDn.*(ec0 - ec1);
    DepsiloncDtau = DfcDtau.*(ec0 - ec1);
    Vc1 = [ec + rho.*DepsiloncDnup, ec + rho.*DepsiloncDndn];
    Vc2 = rho.*DepsiloncDDn;
    Vc3 = rho.*DepsiloncDtau;
end