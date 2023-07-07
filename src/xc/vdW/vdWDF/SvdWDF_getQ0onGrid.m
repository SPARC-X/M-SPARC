function [S] = SvdWDF_getQ0onGrid(S, rho, ecPW, v_cPW)
% @file    SvdWDF_getQ0onGrid.m
% @brief   This file contains the functions computing the energy ratio q0(x)
% @authors Boqin Zhang <bzhang376@gatech.edu>
%          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% Reference:
% Dion, Max, Henrik Rydberg, Elsebeth Schröder, David C. Langreth, and Bengt I. Lundqvist. 
% "Van der Waals density functional for general geometries." 
% Physical review letters 92, no. 24 (2004): 246401.
% Román-Pérez, Guillermo, and José M. Soler. 
% "Efficient implementation of a van der Waals density functional: application to double-wall carbon nanotubes." 
% Physical review letters 103, no. 9 (2009): 096102.
% Thonhauser, T., S. Zuluaga, C. A. Arter, K. Berland, E. Schröder, and P. Hyldgaard. 
% "Spin signature of nonlocal correlation binding in metal-organic frameworks." 
% Physical review letters 115, no. 13 (2015): 136402.
% Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
% ==============================================================================================
    nnr = S.Nx*S.Ny*S.Nz;
    q_cut = S.vdWDF_qmesh(end);
    q_min = S.vdWDF_qmesh(1);
%     qnum = size(S.vdWDF_qmesh, 1);
    S.vdWDF_q0 = ones(nnr,1)*q_cut;
    S.vdWDF_Dq0Drho = zeros(nnr, 2); % up and down
    S.vdWDF_Dq0Dgradrho = zeros(nnr, 2); % up and down
    epsr = 1.0E-12;
    
    boolRhoGepsr = rho(:, 1) > epsr;
    boolRhoUpGepsr = (rho(:, 2) > epsr/2.0) & boolRhoGepsr;
    boolRhoDnGepsr = (rho(:, 3) > epsr/2.0) & boolRhoGepsr;
    rhoGepsr = rho(boolRhoGepsr, 1);
    rhoUpGepsr = rho(boolRhoUpGepsr, 2);
    rhoDnGepsr = rho(boolRhoDnGepsr, 3);
    Drho_a1 = S.grad_1*(rho);
    Drho_a2 = S.grad_2*(rho);
    Drho_a3 = S.grad_3*(rho);
    
    directDrhoUp = [Drho_a1(:, 2), Drho_a2(:, 2), Drho_a3(:, 2)];
    DrhoUp_car = S.grad_T'*directDrhoUp';
    S.DrhoUp = DrhoUp_car';
    gradRhoLengthUp = vecnorm(S.DrhoUp, 2, 2);
    
    directDrhoDn = [Drho_a1(:, 3), Drho_a2(:, 3), Drho_a3(:, 3)];
    DrhoDn_car = S.grad_T'*directDrhoDn';
    S.DrhoDn = DrhoDn_car';
    gradRhoLengthDn = vecnorm(S.DrhoDn, 2, 2);
    
    q0x_up = zeros(nnr, 1);
    DqxDrho_up = zeros(nnr, 1);
    q0x_dn = zeros(nnr, 1);
    DqxDrho_dn = zeros(nnr, 1);
    %% q0x_up(dn) is only computed for grids whose rho_up(dn) is more than epsr/2
    fac = 2.0^(-1.0/3.0);
    sUp = gradRhoLengthUp(boolRhoUpGepsr) ./ (2.0*kF(rhoUpGepsr).*rhoUpGepsr);
    qx_up = kF(2.0*rhoUpGepsr) .* Fs(fac*sUp, S.vdWDFFlag);
    [q0xGepsr_up, Dq0xDqxGepsr_up] = saturate_q(qx_up, 4.0*q_cut); % force q to be in an interval
    DqxDrhoGepsr_up = 2.0*Dq0xDqxGepsr_up.*rhoUpGepsr.*dqx_drho(2.0*rhoUpGepsr, fac*sUp, S.vdWDFFlag)...
        + q0xGepsr_up.*rho(boolRhoUpGepsr, 3)./rho(boolRhoUpGepsr, 1);
    
    sDn = gradRhoLengthDn(boolRhoDnGepsr) ./ (2.0*kF(rhoDnGepsr).*rhoDnGepsr);
    qx_dn = kF(2.0*rhoDnGepsr) .* Fs(fac*sDn, S.vdWDFFlag);
    [q0xGepsr_dn, Dq0xDqxGepsr_dn] = saturate_q(qx_dn, 4.0*q_cut); % force q to be in an interval
    DqxDrhoGepsr_dn = 2.0*Dq0xDqxGepsr_dn.*rhoDnGepsr.*dqx_drho(2.0*rhoDnGepsr, fac*sDn, S.vdWDFFlag)...
        + q0xGepsr_dn.*rho(boolRhoDnGepsr, 2)./rho(boolRhoDnGepsr, 1);
    q0x_up(boolRhoUpGepsr) = q0xGepsr_up(:);
    DqxDrho_up(boolRhoUpGepsr) = DqxDrhoGepsr_up(:);
    DqxDrho_up(boolRhoDnGepsr) = DqxDrho_up(boolRhoDnGepsr) - q0xGepsr_dn.*rhoDnGepsr./rho(boolRhoDnGepsr, 1);
    q0x_dn(boolRhoDnGepsr) = q0xGepsr_dn(:);
    DqxDrho_dn(boolRhoDnGepsr) = DqxDrhoGepsr_dn(:);
    DqxDrho_dn(boolRhoUpGepsr) = DqxDrho_dn(boolRhoUpGepsr) - q0xGepsr_up.*rhoUpGepsr./rho(boolRhoUpGepsr, 1);
    %% q is only computed for grids whose rho is larger than epsr
    qxGepsr = (rho(boolRhoGepsr, 2).*q0x_up(boolRhoGepsr) + rho(boolRhoGepsr, 3).*q0x_dn(boolRhoGepsr)) ./ rhoGepsr;
    qcGepsr = -4.0*pi/3.0 * ecPW(boolRhoGepsr);
    qGepsr = qxGepsr + qcGepsr;
    [q0Gepsr, Dq0GepsrDq] = saturate_q(qGepsr, q_cut);
    q0Gepsr(q0Gepsr < q_min) = q_min;
    S.vdWDF_q0(boolRhoGepsr) = q0Gepsr(:);
    Dq0Dq = zeros(nnr, 1);
    Dq0Dq(boolRhoGepsr) = Dq0GepsrDq;
    DqcDrho = -4.0*pi/3.0 * (v_cPW - [ecPW, ecPW]);
    S.vdWDF_Dq0Drho = Dq0Dq.*(DqcDrho + [DqxDrho_up, DqxDrho_dn]); % attention: Dq0Drho at here is (n_up + n_dn)*[d(q0)/d(n_up), d(q0)/d(n_dn)]!
    S.vdWDF_Dq0Dgradrho(boolRhoUpGepsr, 1) = 2.0*Dq0Dq(boolRhoUpGepsr).*Dq0xDqxGepsr_up.*rhoUpGepsr.*kF(2.0*rhoUpGepsr)...
        .*dFs_ds(fac*sUp, S.vdWDFFlag).*ds_dgradrho(2.0*rhoUpGepsr); % attention: Dq0Dgradrho at here is (n_up + n_dn)*d(q0)/d(|grad n_up|)!
    S.vdWDF_Dq0Dgradrho(boolRhoDnGepsr, 2) = 2.0*Dq0Dq(boolRhoDnGepsr).*Dq0xDqxGepsr_dn.*rhoDnGepsr.*kF(2.0*rhoDnGepsr)...
        .*dFs_ds(fac*sDn, S.vdWDFFlag).*ds_dgradrho(2.0*rhoDnGepsr);
end


function [q0, Dq0Dq] = saturate_q(q, q_cutoff)
    nnr = size(q, 1);
    m_cut = 12;
    e_exp = zeros(nnr, 1);
    Dq0Dq = zeros(nnr, 1);
    qDev_qCut = q/q_cutoff;
    for idx = 1:m_cut
        e_exp = e_exp + qDev_qCut.^idx/idx;
        Dq0Dq = Dq0Dq + qDev_qCut.^(idx - 1);
    end
    q0 = q_cutoff*(1.0 - exp(-e_exp));
    Dq0Dq = Dq0Dq.*exp(-e_exp);
end

function kFResult = kF(rho)
    kFResult = (3*pi^2*rho(:)).^(1/3);
end

function FsResult = Fs(s, vdWflag)
    if vdWflag == 1 % vdW-DF1
        Z_ab = -0.8491;
    end
    if vdWflag == 2 %vdW-DF2
        Z_ab = -1.887;
    end
    FsResult = 1.0 - Z_ab*s(:).^2/9.0;
end

function DkFDrho = dkF_drho(rho)
    DkFDrho = (1.0/3.0)*kF(rho)./rho;
end

function DFsDs = dFs_ds(s, vdWflag)
    if vdWflag == 1 % vdW-DF1
        Z_ab = -0.8491;
    end
    if vdWflag == 2 %vdW-DF2
        Z_ab = -1.887;
    end
    DFsDs =  (-2.0/9.0)*s*Z_ab;
end

function DsDrho = ds_drho(rho, s)
    DsDrho = -s.*(dkF_drho(rho)./kF(rho) + 1.0./rho);
end

function DqxDrho = dqx_drho(rho, s, vdWflag)
    DqxDrho = dkF_drho(rho).*Fs(s, vdWflag) + kF(rho).*dFs_ds(s, vdWflag).*ds_drho(rho, s);
end

function DsDgradrho = ds_dgradrho(rho)
    DsDgradrho = 0.5./(kF(rho).*rho);
end
