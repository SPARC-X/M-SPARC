function [S, potential] = vdWDF_uGenerate_Potential(S, inptDrho, vdWDF_Dq0Drho, vdWDF_Dq0Dgradrho, ps, DpDq0s) 
% @file    vdWDF_uGenerate_Potential.m
% @brief   This file contains the functions for getting u in real space,
%          and generating vdW-DF potential
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
	hPrefactor = zeros(nnr, 1);
    % ps and DpDq0s are computed in function vdWDF_splineInterpolation_energy
% 	[ps, DpDq0s] = spline_interpolation_more(S.vdWDF_qmesh, S.vdWDF_q0, S.vdWDF_D2yDx2);
	u = S.vdW_u; % attention: it is NOT wrong! u at here is thetas*Kernel(Phi)
    Drho = inptDrho;
	Dq0Drho = vdWDF_Dq0Drho;
	Dq0Dgradrho = vdWDF_Dq0Dgradrho;
	potential = zeros(nnr, 1);
    qnum = size(S.vdWDF_qmesh, 1);
	for q = 1:qnum
		potential = potential + u(:, q).*(ps(:, q) + DpDq0s(:, q).*Dq0Drho(:));
		hPrefactor = hPrefactor + u(:, q).*DpDq0s(:, q).*Dq0Dgradrho(:);
	end
    gradRhoLength = sqrt(Drho(:, 1).^2 + Drho(:, 2).^2 + Drho(:, 3).^2);
    gradLarger0 = gradRhoLength > 0;
    S.vdW_DpDq0s = DpDq0s; % to be used in calculating stress
%     S.vdW_gradRhoLength = gradRhoLength; % to be used in calculating stress
%     reciVecLabel = S.reciVecLabel;
%     reciVecCoord = S.reciVecCoord;
    for direc = 1:3
        h = hPrefactor.*Drho(:, direc);
        hGradLarger0 = h(gradLarger0)./gradRhoLength(gradLarger0);
        h(gradLarger0) = hGradLarger0(:);
        % hFFT = fftn(reshape(h, S.Nx,S.Ny,S.Nz)); %*(1/nnr)
        % hFFT = hFFT(:);
        % hFFT(S.reciVecLabel(:)) = 1i*S.reciVecCoord(:, direc)*hFFT(S.reciVecLabel(:));
        % hDeriv = ifftn(reshape(hFFT, S.Nx,S.Ny,S.Nz));
        % potential = potential - real(hDeriv(:));
        if S.cell_typ ~= 2
            if direc == 1
                hDirectDeriv(:, 1) = S.grad_1*h;
            elseif direc == 2
                hDirectDeriv(:, 2) = S.grad_2*h;
            else
                hDirectDeriv(:, 3) = S.grad_3*h;
            end
        else
            Dh_1 = S.grad_1*h;
            Dh_2 = S.grad_2*h;
            Dh_3 = S.grad_3*h;
            directDh = [Dh_1, Dh_2, Dh_3];
            Dh_cartesian = S.grad_T'*directDh';
            Dh_cartesian = Dh_cartesian';
            if direc == 1
                hDirectDeriv(:, 1) = Dh_cartesian(:, 1);
            elseif direc == 2
                hDirectDeriv(:, 2) = Dh_cartesian(:, 2);
            else
                hDirectDeriv(:, 3) = Dh_cartesian(:, 3);
            end
        end
    end
    potential = potential - (sum(hDirectDeriv, 2));
% 	S.vdWpotential = potential;
end
