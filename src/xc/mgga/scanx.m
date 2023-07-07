function [ex, v1x, v2x, v3x] = scanx(rho,sigma,tau)
% @file    scanx.m
% @brief   This file contains the functions computing the exchange energy density \epsilon_x and the potential V
% @authors Boqin Zhang <bzhang376@gatech.edu>
%          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% Reference:
% Sun, Jianwei, Adrienn Ruzsinszky, and John P. Perdew. 
% "Strongly constrained and appropriately normed semilocal density functional." 
% Physical review letters 115, no. 3 (2015): 036402.
% Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
% ==============================================================================================
% rho is electron density n, nnr*1 vector 
% !!! attention: rho(rho < S.xc_rhotol) = S.xc_rhotol; process rho before
% calling the function!
% Drho is gradient of electron density n, nnr*3 vectors
% tau is kinetic energy density, nnr*1 vector
    normDrho = sigma.^0.5;
    [s, alpha, DsDn, DsDDn, DalphaDn, DalphaDDn, DalphaDtau] = basicMGGAvariables(rho, normDrho, tau);
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
        + (b1*s.^2 + b2*(1 - alpha).*exp(-b3*(1 - alpha).^2)).^2;
    hx1 = 1 + k1 - k1./(1 + x/k1);
    % interpolate and extrapolate h_x to get F_x
    hx0 = 1.174;
    % switching function f_x
    c1x = 0.667;
    c2x = 0.8;
    dx = 1.24;
    alphaG1 = (alpha > 1);
    alphaE1 = (alpha == 1);
    alphaL1 = (alpha < 1);
    fx = zeros(size(rho, 1), size(rho, 2));
    fx(alphaG1) = -dx*exp(c2x ./ (1 - alpha(alphaG1)));
    fx(alphaE1) = 0.0;
    fx(alphaL1) = exp(-c1x*alpha(alphaL1) ./ (1 - alpha(alphaL1)));
    a1 = 4.9479;
    gx = 1 - exp(-a1*s.^(-0.5));
    Fx = (hx1 + fx.*(hx0 - hx1)).*gx;
    % get epsilon_x
    ex = epsixUnif.*Fx;
    %% solve variation of F_x
    s2 = s.*s;
    term1 = 1 + (b4*s2)/muak.*exp(-abs(b4)*s2/muak);
    term2 = s2.*(b4/muak*exp(-abs(b4)*s2/muak) + b4*s2/muak.*exp(-abs(b4)*s2/muak)*(-abs(b4)/muak));
    term3 = 2*(b1*s2 + b2*(1 - alpha).*exp(-b3*(1 - alpha).^2));
    term4 = b2*(-exp(-b3*(1-alpha).^2) + (1-alpha).*exp(-b3*(1-alpha).^2).*(2*b3*(1-alpha)));
    DxDs = 2*s.*(muak*(term1 + term2) + b1*term3);
    DxDalpha = term3.*term4;
    DxDn = DsDn.*DxDs + DalphaDn.*DxDalpha;
    DxDDn = DsDDn.*DxDs + DalphaDDn.*DxDalpha;
    DxDtau = DalphaDtau.*DxDalpha;
    
    DgxDn = -exp(-a1*s.^(-0.5)).*(a1/2*s.^(-1.5)).*DsDn;
    DgxDDn = -exp(-a1*s.^(-0.5)).*(a1/2*s.^(-1.5)).*DsDDn;
    Dhx1Dx = 1 ./ (1 + x/k1).^2;
    Dhx1Dn = DxDn .* Dhx1Dx;
    Dhx1DDn = DxDDn .* Dhx1Dx;
    Dhx1Dtau = DxDtau .* Dhx1Dx;
    DfxDalpha = zeros(size(rho, 1), size(rho, 2));
    DfxDalpha(alphaG1) = -dx*exp(c2x./(1-alpha(alphaG1))).*(c2x./(1-alpha(alphaG1)).^2);
    DfxDalpha(alphaE1) = 0.0;
    DfxDalpha(alphaL1) = exp(-c1x*alpha(alphaL1)./(1-alpha(alphaL1))).*(-c1x./(1-alpha(alphaL1)).^2);
    DfxDn = DfxDalpha.*DalphaDn;
    DfxDDn = DfxDalpha.*DalphaDDn;
    DfxDtau = DfxDalpha.*DalphaDtau;
    
    DFxDn = (hx1+fx.*(hx0-hx1)).*DgxDn + gx.*(1-fx).*Dhx1Dn + gx.*(hx0-hx1).*DfxDn;
    DFxDDn = (hx1+fx.*(hx0-hx1)).*DgxDDn + gx.*(1-fx).*Dhx1DDn + gx.*(hx0-hx1).*DfxDDn;
    DFxDtau = gx.*(1-fx).*Dhx1Dtau + gx.*(hx0-hx1).*DfxDtau;
    %% solve variant of n*epsilon_x^{unif}*F_x
    DepsixUnifDn = -(3*pi^2)^(1/3)/(4*pi) * rho.^(-2/3);
    v1x = (epsixUnif+rho.*DepsixUnifDn).*Fx + rho.*epsixUnif.*DFxDn;
    v2x = rho.*epsixUnif.*DFxDDn./normDrho; % D(rho*Ex)/D(|grad rho|) / |grad rho|
    v3x = rho.*epsixUnif.*DFxDtau;
end