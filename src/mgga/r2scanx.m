function [ex,vx,v2x,v3x] = r2scanx(rho,sigma,tau)
% @file    r2scanx.m
% @brief   This file contains the functions computing the exchange energy density \epsilon_x and the potential V
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
    [p, alpha, DpDn, DpDDn, DalphaDn, DalphaDDn, DalphaDtau] = ex_basicr2SCANvariables(rho, normDrho, tau);
    [ex, vx, v2x, v3x] = exchanger2SCAN(rho, p, alpha, DpDn, DpDDn, DalphaDn, DalphaDDn, DalphaDtau);
    v2x = v2x ./ normDrho;  % D(rho*Ex)/D(|grad rho|) / |grad rho|
end

function [p, alpha, DpDn, DpDDn, DalphaDn, DalphaDDn, DalphaDtau] ...
    = ex_basicr2SCANvariables(rho, normDrho, tau)
%     normDrho = (Drho(:,1).^2 + Drho(:,2).^2 + Drho(:,3).^2).^0.5;
    %% value
    s = normDrho ./ (2*(3*pi^2)^(1/3)*rho.^(4/3));
    p = s .^2;
    tauw = normDrho.^2 ./ (8*rho);
    tauUnif = 3/10*(3*pi^2)^(2/3) * rho.^(5/3);
    eta = 0.001;
    alpha = (tau - tauw) ./ (tauUnif + eta*tauw); % when there is no spin, ds = 1
    %% variation
    DsDn = -2*normDrho ./ (3*(3*pi^2)^(1/3).*rho.^(7/3));
    DpDn = 2*s .* DsDn;
    DsDDn = 1 ./ (2*(3*pi^2)^(1/3)*rho.^(4/3));
    DpDDn = 2*s .* DsDDn;
    DtauwDn = -normDrho.^2 ./ (8*rho.^2);
    DtauwDDn = normDrho ./ (4*rho);
    DtauUnifDn = (3*pi^2)^(2/3)/2 * rho.^(2/3);
    DalphaDn = (-DtauwDn.*(tauUnif + eta*tauw) - (tau - tauw).*(DtauUnifDn + eta*DtauwDn)) ./ ((tauUnif + eta*tauw).^2);
    DalphaDDn = (-DtauwDDn.*(tauUnif + eta*tauw) - (tau - tauw).*eta.*DtauwDDn) ./ ((tauUnif + eta*tauw).^2);
    DalphaDtau = 1 ./ (tauUnif + eta*tauw);
end

function [ex, vx, v2x, v3x] = exchanger2SCAN(rho, p, alpha, DpDn, DpDDn, DalphaDn, DalphaDDn, DalphaDtau)
    %% solve epsilon_x
    epsixUnif = -3/(4*pi)*(3*pi^2*rho).^(1/3);
    % compose h_x^1
    k0 = 0.174;
    k1 = 0.065;
    muak = 10/81;
    % compose x
    eta = 0.001;
    Ceta = 20/27 + eta*5/3;
    C2 = -0.162742;
    dp2 = 0.361;
    x = (Ceta*C2 * exp(-p.^2 / dp2^4) + muak) .*p;
    hx0 = 1 + k0;
    hx1 = 1 + k1 - k1./(1 + x/k1);
    % interpolate and extrapolate h_x to get F_x
    % switching function f_x
    c1x = 0.667;
    c2x = 0.8;
    dx = 1.24;
    alphaG25 = (alpha > 2.5);
    alpha0To25 = ((alpha >= 0.0) & (alpha <= 2.5));
    alphaL0 = (alpha < 0.0);
    fx = zeros(size(rho, 1), size(rho, 2));
    fx(alphaG25) = -dx*exp(c2x ./ (1 - alpha(alphaG25)));
    fx(alpha0To25) = 1.0 + (-0.667)*alpha(alpha0To25) + (-0.4445555)*alpha(alpha0To25).^2 + (-0.663086601049)*alpha(alpha0To25).^3 ...
        + 1.451297044490*alpha(alpha0To25).^4 + (-0.887998041597)*alpha(alpha0To25).^5 + 0.234528941479*alpha(alpha0To25).^6 ...
        + (-0.023185843322)*alpha(alpha0To25).^7;
    fx(alphaL0) = exp(-c1x*alpha(alphaL0) ./ (1 - alpha(alphaL0)));
    a1 = 4.9479;
    gx = 1 - exp(-a1*p.^(-0.25));
    Fx = (hx1 + fx.*(hx0 - hx1)).*gx;
    % get epsilon_x
    ex = epsixUnif.*Fx;
    %% solve variation of F_x
    DxDp = (Ceta*C2 * exp(-p.^2 / dp2^4) + muak) + Ceta*C2*exp(-p.^2 / dp2^4).*( -2*p/ dp2^4).*p;
    DxDn = DxDp.*DpDn;
    DxDDn = DpDDn.*DxDp;
    
    DgxDn = -exp(-a1*p.^(-0.25)).*(a1/4*p.^(-1.25)).*DpDn;
    DgxDDn = -exp(-a1*p.^(-0.25)).*(a1/4*p.^(-1.25)).*DpDDn;
    Dhx1Dx = 1 ./ (1 + x/k1).^2;
    Dhx1Dn = DxDn .* Dhx1Dx;
    Dhx1DDn = DxDDn .* Dhx1Dx;
    DfxDalpha = zeros(size(rho, 1), size(rho, 2));
    DfxDalpha(alphaG25) = -dx*exp(c2x./(1-alpha(alphaG25))).*(c2x./(1-alpha(alphaG25)).^2);
    DfxDalpha(alpha0To25) = (-0.667) + (-0.4445555)*alpha(alpha0To25)*2 + (-0.663086601049)*alpha(alpha0To25).^2*3 ...
        + 1.451297044490*alpha(alpha0To25).^3*4 + (-0.887998041597)*alpha(alpha0To25).^4*5 + 0.234528941479*alpha(alpha0To25).^5*6 ...
        + (-0.023185843322)*alpha(alpha0To25).^6*7;
    DfxDalpha(alphaL0) = exp(-c1x*alpha(alphaL0)./(1-alpha(alphaL0))).*(-c1x./(1-alpha(alphaL0)).^2);
    DfxDn = DfxDalpha.*DalphaDn;
    DfxDDn = DfxDalpha.*DalphaDDn;
    DfxDtau = DfxDalpha.*DalphaDtau;
    
    DFxDn = (hx1+fx.*(hx0-hx1)).*DgxDn + gx.*(1-fx).*Dhx1Dn + gx.*(hx0-hx1).*DfxDn;
    DFxDDn = (hx1+fx.*(hx0-hx1)).*DgxDDn + gx.*(1-fx).*Dhx1DDn + gx.*(hx0-hx1).*DfxDDn;
    DFxDtau = gx.*(hx0-hx1).*DfxDtau;
    %% solve variant of n*epsilon_x^{unif}*F_x
    DepsixUnifDn = -(3*pi^2)^(1/3)/(4*pi) * rho.^(-2/3);
    vx = (epsixUnif+rho.*DepsixUnifDn).*Fx + rho.*epsixUnif.*DFxDn;
    v2x = rho.*epsixUnif.*DFxDDn; % to compare it with Libxc, Vx2/2/|grad rho|
    v3x = rho.*epsixUnif.*DFxDtau;
end