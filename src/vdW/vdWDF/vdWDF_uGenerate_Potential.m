function S = vdWDF_uGenerate_Potential(S)
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
	qnum = size(S.vdWDF_qmesh, 1);
    for q = 1:qnum
        uReci3D = reshape(S.vdW_u(:, q), S.Nx, S.Ny, S.Nz);
        u3D = ifftn(uReci3D)*nnr; % the position of const 1/(Nx*Ny*Nz) is in ifftn in MATLAB 
        %but the const is in fwfft(cfft3d R->G) in FORTRAN
        S.vdW_u(:, q) = u3D(:);
    end
	hPrefactor = zeros(nnr, 1);
	[ps, DpDq0s] = spline_interpolation_more(S.vdWDF_qmesh, S.vdWDF_q0, S.vdWDF_D2yDx2);
	u = S.vdW_u; % attention: it is NOT wrong! u at here is thetas*Kernel(Phi)
	Dq0Drho = S.vdWDF_Dq0Drho;
	Dq0Dgradrho = S.vdWDF_Dq0Dgradrho;
	potential = zeros(nnr, 1);
	for q = 1:qnum
		potential = potential + u(:, q).*(ps(:, q) + DpDq0s(:, q).*Dq0Drho(:));
		hPrefactor = hPrefactor + u(:, q).*DpDq0s(:, q).*Dq0Dgradrho(:);
	end
    gradRhoLength = sqrt(S.Drho(:, 1).^2 + S.Drho(:, 2).^2 + S.Drho(:, 3).^2);
    gradLarger0 = gradRhoLength > 0;
    S.vdW_DpDq0s = DpDq0s; % to be used in calculating stress
    S.vdW_gradRhoLength = gradRhoLength; % to be used in calculating stress
%     reciVecLabel = S.reciVecLabel;
%     reciVecCoord = S.reciVecCoord;
    for direc = 1:3
        h = hPrefactor.*S.Drho(:, direc);
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
	S.vdWpotential = potential;
end

function [ps, DpDq0s] = spline_interpolation_more(qmesh, q0, D2yDx2)
    qnum = size(qmesh, 1);
    gridNum = size(q0, 1);
    ps = zeros(gridNum, qnum);
    DpDq0s = zeros(gridNum, qnum);
    lowerBound = ones(gridNum, 1);
    upperBound = ones(gridNum, 1)*qnum;
    while (max(upperBound - lowerBound) > 1)
        idx = floor((upperBound + lowerBound)/2);
        boolLargerq0idx = q0 > qmesh(idx(:));
        lowerBound(boolLargerq0idx) = idx(boolLargerq0idx);
        upperBound(~boolLargerq0idx) = idx(~boolLargerq0idx);
    end
    dx = qmesh(upperBound) - qmesh(lowerBound);
    a = (qmesh(upperBound) - q0)./dx;
    b = (q0 - qmesh(lowerBound))./dx;
    c = (a.^3 - a).*(dx.^2)/6.0;
    d = (b.^3 - b).*(dx.^2)/6.0;
    e = (3.0*a.^2 - 1.0).*dx/6.0;
    f = (3.0*b.^2 - 1.0).*dx/6.0;
    for q = 1:qnum
        y = zeros(qnum, 1);
        y(q) = 1;
        for i=1:gridNum
            ps(i, q) = a(i)*y(lowerBound(i)) + b(i)*y(upperBound(i));
            ps(i, q) = ps(i, q) + c(i)*D2yDx2(q, lowerBound(i));
            ps(i, q) = ps(i, q) + d(i)*D2yDx2(q, upperBound(i));
            DpDq0s(i, q) = (y(upperBound(i)) - y(lowerBound(i)))/dx(i) - e(i)*D2yDx2(q,lowerBound(i)) + f(i)*D2yDx2(q,upperBound(i));
        end
    end
end