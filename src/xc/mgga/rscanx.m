function [ex, vx, v2x, v3x] = rscanx(rho,sigma,tau)
% @file    rscanx.m
% @brief   This file contains the functions computing the exchange energy density \epsilon_x and the potential V
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
    [ex, vx, v2x, v3x] = exchangerSCAN(rho, s, alphaP, DsDn, DsDDn, DalphaPDn, DalphaPDDn, DalphaPDtau);
    v2x = v2x ./ normDrho;  % D(rho*Ex)/D(|grad rho|) / |grad rho|
end

function [ex, vx, v2x, v3x] = exchangerSCAN(rho, s, alphaP, DsDn, DsDDn, DalphaPDn, DalphaPDDn, DalphaPDtau)
    %% solve epsilon_x
    epsixUnif = -3/(4*pi)*(3*pi^2*rho).^(1/3);
    % compose h_x^1
    k1 = 0.065;
    muak = 10/81;
    b2 = sqrt(5913/405000);
    b1 = 511/13500/(2*b2);
    b3 = 0.5;
    b4 = muak*muak/k1 - 1606/18225 - b1*b1;
    % compose x
    x = muak*s.^2 .* (1 + b4*s.^2/muak .* exp(-abs(b4)*s.^2/muak))...
        + (b1*s.^2 + b2*(1 - alphaP).*exp(-b3*(1 - alphaP).^2)).^2;
    hx1 = 1 + k1 - k1./(1 + x/k1);
    % interpolate and extrapolate h_x to get F_x
    hx0 = 1.174;
    % switching function f_x
    c1x = 0.667;
    c2x = 0.8;
    dx = 1.24;
    alphaG25 = (alphaP > 2.5);
    alpha0To25 = ((alphaP >= 0.0) & (alphaP <= 2.5));
    alphaL0 = (alphaP < 0.0);
    fx = zeros(size(rho, 1), size(rho, 2));
    fx(alphaG25) = -dx*exp(c2x ./ (1 - alphaP(alphaG25)));
    fx(alpha0To25) = 1.0 + (-0.667)*alphaP(alpha0To25) + (-0.4445555)*alphaP(alpha0To25).^2 + (-0.663086601049)*alphaP(alpha0To25).^3 ...
        + 1.451297044490*alphaP(alpha0To25).^4 + (-0.887998041597)*alphaP(alpha0To25).^5 + 0.234528941479*alphaP(alpha0To25).^6 ...
        + (-0.023185843322)*alphaP(alpha0To25).^7;
    fx(alphaL0) = exp(-c1x*alphaP(alphaL0) ./ (1 - alphaP(alphaL0)));
    a1 = 4.9479;
    gx = 1 - exp(-a1*s.^(-0.5));
    Fx = (hx1 + fx.*(hx0 - hx1)).*gx;
    % get epsilon_x
    ex = epsixUnif.*Fx;
    %% solve variation of F_x
    s2 = s.*s;
    term1 = 1 + (b4*s2)/muak.*exp(-abs(b4)*s2/muak);
    term2 = s2.*(b4/muak*exp(-abs(b4)*s2/muak) + b4*s2/muak.*exp(-abs(b4)*s2/muak)*(-abs(b4)/muak));
    term3 = 2*(b1*s2 + b2*(1 - alphaP).*exp(-b3*(1 - alphaP).^2));
    term4 = b2*(-exp(-b3*(1-alphaP).^2) + (1-alphaP).*exp(-b3*(1-alphaP).^2).*(2*b3*(1-alphaP)));
    DxDs = 2*s.*(muak*(term1 + term2) + b1*term3);
    DxDalpha = term3.*term4;
    DxDn = DsDn.*DxDs + DalphaPDn.*DxDalpha;
    DxDDn = DsDDn.*DxDs + DalphaPDDn.*DxDalpha;
    DxDtau = DalphaPDtau.*DxDalpha;
    
    DgxDn = -exp(-a1*s.^(-0.5)).*(a1/2*s.^(-1.5)).*DsDn;
    DgxDDn = -exp(-a1*s.^(-0.5)).*(a1/2*s.^(-1.5)).*DsDDn;
    Dhx1Dx = 1 ./ (1 + x/k1).^2;
    Dhx1Dn = DxDn .* Dhx1Dx;
    Dhx1DDn = DxDDn .* Dhx1Dx;
    Dhx1Dtau = DxDtau .* Dhx1Dx;
    DfxDalpha = zeros(size(rho, 1), size(rho, 2));
    DfxDalpha(alphaG25) = -dx*exp(c2x./(1-alphaP(alphaG25))).*(c2x./(1-alphaP(alphaG25)).^2);
    DfxDalpha(alpha0To25) = (-0.667) + (-0.4445555)*alphaP(alpha0To25)*2 + (-0.663086601049)*alphaP(alpha0To25).^2*3 ...
        + 1.451297044490*alphaP(alpha0To25).^3*4 + (-0.887998041597)*alphaP(alpha0To25).^4*5 + 0.234528941479*alphaP(alpha0To25).^5*6 ...
        + (-0.023185843322)*alphaP(alpha0To25).^6*7;
    DfxDalpha(alphaL0) = exp(-c1x*alphaP(alphaL0)./(1-alphaP(alphaL0))).*(-c1x./(1-alphaP(alphaL0)).^2);
    DfxDn = DfxDalpha.*DalphaPDn;
    DfxDDn = DfxDalpha.*DalphaPDDn;
    DfxDtau = DfxDalpha.*DalphaPDtau;
    
    DFxDn = (hx1+fx.*(hx0-hx1)).*DgxDn + gx.*(1-fx).*Dhx1Dn + gx.*(hx0-hx1).*DfxDn;
    DFxDDn = (hx1+fx.*(hx0-hx1)).*DgxDDn + gx.*(1-fx).*Dhx1DDn + gx.*(hx0-hx1).*DfxDDn;
    DFxDtau = gx.*(1-fx).*Dhx1Dtau + gx.*(hx0-hx1).*DfxDtau;
    %% solve variant of n*epsilon_x^{unif}*F_x
    DepsixUnifDn = -(3*pi^2)^(1/3)/(4*pi) * rho.^(-2/3);
    vx = (epsixUnif+rho.*DepsixUnifDn).*Fx + rho.*epsixUnif.*DFxDn;
    v2x = rho.*epsixUnif.*DFxDDn; % D(rho*Ex)/D(|grad rho|)
    v3x = rho.*epsixUnif.*DFxDtau;
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
    
