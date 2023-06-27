function [ex,vx,v2x,v3x] = scanx_spin(rho,sigma,tau)
% @file    scanx_spin.m
% @brief   This file contains the functions computing the spin-polarized exchange energy density \epsilon_x and the potential V
% @authors Boqin Zhang <bzhang376@gatech.edu>
%          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% Reference:
% Sun, Jianwei, Adrienn Ruzsinszky, and John P. Perdew. 
% "Strongly constrained and appropriately normed semilocal density functional." 
% Physical review letters 115, no. 3 (2015): 036402.
% Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
% ==============================================================================================
 % all three input variables have 3 cols
 % rho col1: n; col2: n up; col3: n dn
 % sigma col1: grad of n; col2: grad of n up; col3: grad of n dn
 % ex has only 1 column: epsilon_x = [\epsilon_x(2n_up)*2n_up + \epsilon_x(2n_dn)*2n_dn] / 2 / (n_up + n_dn)
 % vx has 2 cols: [D(\epsilon_x(2n_up)*2n_up)/D(2n_up), D(\epsilon_x(2n_dn)*2n_dn)/D(2n_dn)]
 % v2x has 2 cols: [D(\epsilon_x(2n_up)*2n_up)/D(|grad 2n_up|) / |grad 2n_up|*2, D(\epsilon_x(2n_dn)*2n_dn)/D(|grad 2n_dn|) / |grad 2n_dn|*2]
 %                =[D(\epsilon_x(2n_up)*n_up)/D(|grad n_up|) / |grad n_up|, D(\epsilon_x(2n_dn)*n_dn)/D(|grad n_dn|) / |grad n_dn|]
 % v3x has 2 cols: [D(\epsilon_x(2n_up)*2n_up)/D(2tau_up), D(\epsilon_x(2n_dn)*2n_dn)/D(2tau_dn)]
    [ex_up, Vx1_up, Vx2_up, Vx3_up] = scanx(2*rho(:, 2), 4*sigma(:, 2), 2*tau(:, 2));
    [ex_dn, Vx1_dn, Vx2_dn, Vx3_dn] = scanx(2*rho(:, 3), 4*sigma(:, 3), 2*tau(:, 3));
    ex = (ex_up.*rho(:, 2) + ex_dn.*rho(:, 3)) ./ rho(:, 1);
    vx = [Vx1_up, Vx1_dn];
    v2x = [Vx2_up, Vx2_dn] * 2.0;
    v3x = [Vx3_up, Vx3_dn];
end