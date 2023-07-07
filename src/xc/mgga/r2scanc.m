function [ec,vc,v2c,v3c] = r2scanc(rho,sigma,tau)
% @file    r2scanx.m
% @brief   This file contains the functions computing the correlation energy density \epsilon_x and the potential V
% @authors Boqin Zhang <bzhang376@gatech.edu>
%          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% Reference:
% Furness, James W., Aaron D. Kaplan, Jinliang Ning, John P. Perdew, and Jianwei Sun. 
% "Accurate and numerically efficient r2SCAN meta-generalized gradient approximation." 
% The journal of physical chemistry letters 11, no. 19 (2020): 8208-8215.
% Physical review letters 115, no. 3 (2015): 036402.
% Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
% ==============================================================================================
    normDrho = sigma.^0.5;
    [s, p, alpha, DsDn, DsDDn, DpDn, DpDDn, DalphaDn, DalphaDDn, DalphaDtau] = co_basicr2SCANvariables(rho, normDrho, tau);
    [ec, vc, v2c, v3c] = correlationr2SCAN(rho, s, p, alpha, DsDn, DsDDn, DpDn, DpDDn, DalphaDn, DalphaDDn, DalphaDtau);
    v2c = v2c ./ normDrho; % D(rho*Ec)/D(|grad rho|) / |grad rho|
end

function [s, p, alpha, DsDn, DsDDn, DpDn, DpDDn, DalphaDn, DalphaDDn, DalphaDtau] ...
    = co_basicr2SCANvariables(rho, normDrho, tau)
    %% value
    s = normDrho ./ (2*(3*pi^2)^(1/3)*rho(:, 1).^(4/3));
    p = s .^2;

    zeta = 0.0;
    tauw = normDrho.^2 ./ (8*rho(:, 1));

    ds = 1.0; % when there is no spin, ds = 1
    tauUnif = 3/10*(3*pi^2)^(2/3) * rho(:, 1).^(5/3);
    eta = 0.001;
    alpha = (tau - tauw) ./ (tauUnif + eta*tauw); % when there is no spin, ds = 1
    %% variation
    DsDn = -2*normDrho ./ (3*(3*pi^2)^(1/3).*rho(:, 1).^(7/3));
    DpDn = 2*s .* DsDn;
    DsDDn = 1 ./ (2*(3*pi^2)^(1/3)*rho(:, 1).^(4/3));
    DpDDn = 2*s .* DsDDn;
    DtauwDn = -normDrho.^2 ./ (8*rho(:, 1).^2);
    DtauwDDn = normDrho ./ (4*rho(:, 1));

    DtauUnifDn = (3*pi^2)^(2/3)/2 * rho(:, 1).^(2/3);

    DalphaDn = (-DtauwDn.*(tauUnif + eta*tauw) - (tau - tauw).*(DtauUnifDn + eta*DtauwDn)) ./ ((tauUnif + eta*tauw).^2);
    DalphaDDn = (-DtauwDDn.*(tauUnif + eta*tauw) - (tau - tauw).*eta.*DtauwDDn) ./ ((tauUnif + eta*tauw).^2);
    DalphaDtau = 1 ./ (tauUnif + eta*tauw);
end

function [ec, vc, v2c, v3c] = correlationr2SCAN(rho, s, p, alpha, DsDn, DsDDn, DpDn, DpDDn, DalphaDn, DalphaDDn, DalphaDtau)
    % all input variables of this function only have one column. 
    %     phi = ((1 + zeta).^(2/3) + (1 - zeta).^(2/3)) / 2;
    phi = 1.0;
    
    rs = (0.75./(pi*rho)).^(1/3) ;
    %% compute epsilon_c^0 when alpha \approx 0
    b1c = 0.0285764;
    b2c = 0.0889;
    b3c = 0.125541;
    ecLDA0 = -b1c ./ (1 + b2c*(rs.^0.5) + b3c*rs);
    dx = 1;
    cx0 = -3/(4*pi)*(9*pi/4)^(1/3); % unused variable
    Gc = 1.0;
    w0 = exp(-ecLDA0/b1c) - 1;
    betaConst = 0.06672455060314922;
    betaRsInf = betaConst*0.1/0.1778; % the value is estimated
    f0 = -0.9;
    chiInf = (3*pi^2/16)^(2/3) * (betaRsInf*1/(cx0 - f0)); % the value is estimated, % should be 0.128026 if zeta = 0
    gInf0s = (1 + 4*chiInf.*(s.^2)).^(-0.25);
    H0 = b1c*log(1 + w0.*(1 - gInf0s));
    ec0 = ecLDA0 + H0;
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
    
    f_zeta = 0.0;
    fp_zeta = 0.0;
    zeta4 = 0.0;
    
    ec_lsda1 = ecrs0;
    declsda1_drs = decrs0_drs;
    
    r = 0.031091;
    w1 = exp(-ec_lsda1 ./ r) - 1;
    t = (3*pi^2/16)^(1/3) * s./sqr_rs;

    y = beta./(r*w1) .* (t.^2); % new variable, y = A*t*t
    deltafc2 = 1*(-0.64) + 2*(-0.4352) + 3*(-1.535685604549) + 4*3.061560252175 + 5*(-1.915710236206) + 6*0.516884468372 + 7*(-0.051848879792); 
    % new variable, a constant
    ds = 1.0;
    eta = 0.001;
    dp2 = 0.361;
    ec_lsda0 = ecLDA0; % new variable, decLDA0drs is in the formula of DecLDA0Dn
    declsda0_drs = b1c*(0.5*b2c*rs.^(-0.5) + b3c)./((1 + b2c*rs.^0.5 + b3c*rs).^2); % new variable
      
    deltayPart1 = deltafc2 ./ (27*r*w1);
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
    DrsDn = -4*pi/9*(4*pi/3*rho).^(-4/3);
    DgInf0sDs = -0.25*(1 + 4*chiInf*s.*s).^(-1.25) .* (4*chiInf*2*s);
    DgInf0sDn = DgInf0sDs .* DsDn;
    DgInf0sDDn = DgInf0sDs .* DsDDn;
    DecLDA0Dn = b1c*(0.5*b2c*rs.^(-0.5) + b3c)./((1 + b2c*rs.^0.5 + b3c*rs).^2).*DrsDn;
    Dw0Dn = (w0 + 1) .* (-DecLDA0Dn/b1c);
    DH0Dn = b1c*(Dw0Dn.*(1 - gInf0s) - w0.*DgInf0sDn) ./ (1 + w0.*(1 - gInf0s));
    DH0DDn = b1c* (-w0.*DgInf0sDDn) ./ (1 + w0.*(1 - gInf0s));
    Dec0Dn = DecLDA0Dn + DH0Dn;
    Dec0DDn = DH0DDn;
    %% compute variation of epsilon_c^1
    dec_lsda1_drs = decrs0_drs;
    
    dec0log_drs = -ec0_q1p .* ec0_den;
    dec0q0_drs = -2.0 * AAec0 * alpha1ec0;
    dec0q1p_drs = AAec0 * ((-0.5)*beta1ec0*rsm1_2./rs + 3.0*0.5*beta3ec0*rsm1_2 + 4.0*beta4ec0);
    dec0den_drs = -(2*ec0_q1.*ec0_q1p + ec0_q1p) ./ (ec0_q1.*ec0_q1 + ec0_q1).^2;
    d2ecrs0_drs2 = -2.0 * AAec0 * alpha1ec0 * dec0log_drs - (dec0q0_drs .* ec0_q1p .* ec0_den + ...
        ec0_q0 .* dec0q1p_drs .* ec0_den + ec0_q0 .* ec0_q1p .* dec0den_drs);
        
    d2eclsda1_drs2 = d2ecrs0_drs2;

    Ddeclsda1_drsDn = d2eclsda1_drs2.*DrsDn;
    % new variable, derives from decrs0_drs, dmacrs_drs, decrs1_drs formula above
    
    Dec_lsda1Dn = (- rs * 1/3 .* dec_lsda1_drs) ./ rho;
    % 	Dec_lsda1Dnup = (- rs * 1/3 .* dec_lsda1_drs - zeta .* dec_lsda1_dzeta + dec_lsda1_dzeta) ./ rho(:, 1);
    % 	Dec_lsda1Dndn = (- rs * 1/3 .* dec_lsda1_drs - zeta .* dec_lsda1_dzeta - dec_lsda1_dzeta) ./ rho(:, 1); % from LDA_PW(S)
    %     fprintf("ecLSDA %10.8f, Vc1LSDA %10.8f %10.8f\n", ec_lsda1, Dec_lsda1Dnup, Dec_lsda1Dndn);
    DbetaDn = 0.066725*(0.1*(1+0.1778*rs) - 0.1778*(1+0.1*rs)) ./ ((1+0.1778*rs).^2) .* DrsDn;
    %     DphiDn = 0.5*(2/3*(1+zeta).^(-1/3) - 2/3*(1-zeta).^(-1/3)).*DzetaDn; % no spin, should be 0
    DtDn = (3*pi^2/16)^(1/3)*(sqr_rs.*DsDn - s.*(DrsDn./(2*sqr_rs))) ./ rs;
    DtDDn = t.*DsDDn./s;

    Dw1Dn = (w1 + 1) .* (-(r*Dec_lsda1Dn) ./ ((r).^2));

    DyDn = (w1.*DbetaDn - beta.*Dw1Dn) ./ (r*w1.^2) .* (t.^2) + beta./(r*w1) .* (2*t).*DtDn; % new variable
    DyDDn = beta./(r*w1) .* (2*t).*DtDDn; % new variable
        
    Declsda0Dn = declsda0_drs.*DrsDn;
    d2eclsda0_drs2 = b1c*((0.5*b2c*(-0.5)*rs.^(-1.5)).*((1 + b2c*rs.^0.5 + b3c*rs).^2) - ...
        (0.5*b2c*rs.^(-0.5) + b3c).*2.*(1 + b2c*rs.^0.5 + b3c*rs).*(0.5*b2c*rs.^(-0.5) + b3c)) ...
        ./ (((1 + b2c*rs.^0.5 + b3c*rs).^2).^2);
    Ddeclsda0_drsDn = d2eclsda0_drs2.*DrsDn;
    
    d_deltayPart1_dn =  0.0 + 0.0 + deltafc2 ./ (27*r) .* (-1).*w1.^(-2) .* Dw1Dn;
    
    d_deltayPart2_dn = 20*(declsda0_drs - declsda1_drs).*DrsDn ...
        + 20*rs.*(Ddeclsda0_drsDn - Ddeclsda1_drsDn) ...
        - 45*eta*(Declsda0Dn - Dec_lsda1Dn);

    d_deltayPart3_dp = exp(-p.^2/dp2^4) + p.*exp(-p.^2/dp2^4).*(-2*p/dp2^4);
        
    DdeltayDn = d_deltayPart1_dn.*deltayPart2.*deltayPart3 ...
        + deltayPart1.*d_deltayPart2_dn.*deltayPart3 ...
        + deltayPart1.*deltayPart2.*d_deltayPart3_dp.*DpDn; % new variable
    DdeltayDDn = deltayPart1.*deltayPart2.*d_deltayPart3_dp.*DpDDn;
    
    DgDn = -0.25*(1 + 4*(y - deltay)).^(-1.25) .* (4*(DyDn - DdeltayDn)); % new formula
    DgDDn = -0.25*(1 + 4*(y - deltay)).^(-1.25) .* (4*(DyDDn - DdeltayDDn));

    DH1Dn = r*((Dw1Dn.*(1 - g) - w1.*DgDn) ./ (1 + w1.*(1 - g)));
    DH1DDn = r* ((-w1.*DgDDn) ./ (1 + w1.*(1 - g)));
    
    Dec1Dn = Dec_lsda1Dn + DH1Dn;
    Dec1DDn = DH1DDn;
    %% variant of fc and ec
    DfcDalpha = zeros(size(rho, 1), size(rho, 2));
    DfcDalpha(alphaG25) = fc(alphaG25).*(c2c./(1 - alpha(alphaG25)).^2);
    DfcDalpha(alpha0To25) = (-0.64) + (-0.4352)*alpha(alpha0To25)*2 + (-1.535685604549)*alpha(alpha0To25).^2*3 ...
        + 3.061560252175*alpha(alpha0To25).^3*4 + (-1.915710236206)*alpha(alpha0To25).^4*5 + 0.516884468372*alpha(alpha0To25).^5*6 ...
        + (-0.051848879792)*alpha(alpha0To25).^6*7; % new formula for the interval;
    DfcDalpha(alphaL0) = fc(alphaL0).*(-c1c./(1 - alpha(alphaL0)).^2);
    DfcDn = DfcDalpha.*DalphaDn;
    DfcDDn = DfcDalpha.*DalphaDDn;
    DfcDtau = DfcDalpha.*DalphaDtau;
    DepsiloncDn = Dec1Dn + fc.*(Dec0Dn -Dec1Dn) + DfcDn.*(ec0 - ec1);
    DepsiloncDDn = Dec1DDn + fc.*(Dec0DDn -Dec1DDn) + DfcDDn.*(ec0 - ec1);
    DepsiloncDtau = DfcDtau.*(ec0 - ec1);
    vc = ec + rho.*DepsiloncDn;
    v2c = rho.*DepsiloncDDn;
    v3c = rho.*DepsiloncDtau;
end