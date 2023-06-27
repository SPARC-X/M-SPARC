function [S] = vdWDF_getQ0onGrid(S, rho, ecPW, v_cPW)
% @file    vdWDF_getQ0onGrid.m
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
% Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
% ==============================================================================================
    nnr = S.Nx*S.Ny*S.Nz;
    q_cut = S.vdWDF_qmesh(end);
    q_min = S.vdWDF_qmesh(1);
%     qnum = size(S.vdWDF_qmesh, 1);
    S.vdWDF_q0 = ones(nnr,1)*q_cut;
    S.vdWDF_Dq0Drho = zeros(nnr, 1);
    S.vdWDF_Dq0Dgradrho = zeros(nnr, 1);
    epsr = 1.0E-12;
    boolRhoGepsr = rho > epsr;
    rhoGepsr = rho(boolRhoGepsr);
    Drho_a1 = S.grad_1*(rho);
    Drho_a2 = S.grad_2*(rho);
    Drho_a3 = S.grad_3*(rho);
    directDrho = [Drho_a1, Drho_a2, Drho_a3];
    Drho_car = S.grad_T'*directDrho';
    Drho_x = Drho_car(1, :)';
    Drho_y = Drho_car(2, :)';
    Drho_z = Drho_car(3, :)';
    
    S.Drho = zeros(nnr, 3);
    S.Drho(:, 1) = Drho_x; S.Drho(:, 2) = Drho_y; S.Drho(:, 3) = Drho_z; 
    Drho_xGepsr = Drho_x(boolRhoGepsr);
    Drho_yGepsr = Drho_y(boolRhoGepsr);
    Drho_zGepsr = Drho_z(boolRhoGepsr);
    %% calculations below are for grids whose rho is more than epsr
    r_s = (3./(4*pi*rhoGepsr)).^(1/3);
    kFResult = kF(rhoGepsr);
    s = sqrt(Drho_xGepsr.^2 + Drho_yGepsr.^2 + Drho_zGepsr.^2)./(2*kFResult.*rhoGepsr);
    % [ecLDA_PW, Dq0Drho_p] = pw(r_s); % Dq0Drho not finish in this step
    FsResult = Fs(s, S.vdWDFFlag);
    q_p = -4.0*pi/3.0*ecPW + kFResult.*FsResult; % energy ratio on every point
    [q0_p, Dq0Dq_p] = saturate_q(q_p, q_cut); % force q to be in an interval
    q0_p(q0_p < q_min) = q_min;
    DqxDrho_p = dqx_drho(rhoGepsr, s, S.vdWDFFlag);
    Dq0Drho_p = Dq0Dq_p.*rhoGepsr.*(-4.0*pi/3.0*(v_cPW - ecPW)./rhoGepsr(:) + DqxDrho_p);
    Dq0Dgradrho_p = Dq0Dq_p.*rhoGepsr.*kFResult.*dFs_ds(s, S.vdWDFFlag).*ds_dgradrho(rhoGepsr);
    %% the real result is the results above adding zeros(for grids whose rho is less than epsr)
    S.vdWDF_q0(boolRhoGepsr) = q0_p(:);
    S.vdWDF_Dq0Drho(boolRhoGepsr) = Dq0Drho_p(:);
    S.vdWDF_Dq0Dgradrho(boolRhoGepsr) = Dq0Dgradrho_p(:);
    clear Drho_xGepsr Drho_yGepsr Drho_zGepsr r_s kFResult s FsResult...
        q0_p DqxDrho_p Dq0Drho_p Dq0Dgradrho_p;
%     pack;
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
