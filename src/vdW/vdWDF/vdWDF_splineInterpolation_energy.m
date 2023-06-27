function [S, ps, DpDq0s] = vdWDF_splineInterpolation_energy(S,rho)
% @file    vdWDF_splineInterpolation_energy.m
% @brief   This file contains the functions for interpolating q0 to get
%          energy ratio cpmponent of every model energy ratio ps(q_1)
%          generating theta=ps*rho, and integrating for getting vdW-DF energy,
%          at last generating u(q2)=\sum_{q1}(kernel(q2,q1)*theta(q1)) in reciprocal space
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
    S.vdWenergy = 0.0;
%% the index of reciprocal lattice mesh grids
    reciLatticeVecs = ((S.lat_uvec.*repmat([S.L1, S.L2, S.L3]', 1, 3)) \ (2*pi * eye(3)))';
    omega = det(S.lat_uvec.*repmat([S.L1, S.L2, S.L3]', 1, 3)); % volume of cell
    reciXlabel2 = [0:floor(S.Nx/2),(floor(S.Nx/2) - S.Nx + 1):-1]; % the time of reciprocal lattice vectors of corresponding entries after FFT
    reciYlabel2 = [0:floor(S.Ny/2),(floor(S.Ny/2) - S.Ny + 1):-1];
    reciZlabel2 = [0:floor(S.Nz/2),(floor(S.Nz/2) - S.Nz + 1):-1];
    [reciXlabel3D2, reciYlabel3D2, reciZlabel3D2] = meshgrid(reciXlabel2, reciYlabel2, reciZlabel2); % the sequence of meshgrid is weird: y-x-z
    reciXlabel3D2 = permute(reciXlabel3D2,[2 1 3]);
    reciYlabel3D2 = permute(reciYlabel3D2,[2 1 3]);
    reciZlabel3D2 = permute(reciZlabel3D2,[2 1 3]);
    ngm = S.Nx*S.Ny*S.Nz;
    reciXcoord = reciXlabel3D2*reciLatticeVecs(1, 1) + reciYlabel3D2*reciLatticeVecs(2, 1) + reciZlabel3D2*reciLatticeVecs(3, 1);
    reciYcoord = reciXlabel3D2*reciLatticeVecs(1, 2) + reciYlabel3D2*reciLatticeVecs(2, 2) + reciZlabel3D2*reciLatticeVecs(3, 2);
    reciZcoord = reciXlabel3D2*reciLatticeVecs(1, 3) + reciYlabel3D2*reciLatticeVecs(2, 3) + reciZlabel3D2*reciLatticeVecs(3, 3);
    vecLength3D = sqrt(reciXcoord.*reciXcoord + reciYcoord.*reciYcoord + reciZcoord.*reciZcoord);
    S.reciVecCoord = zeros(ngm, 3);
    S.reciVecCoord(:, 1) = reciXcoord(:); 
    S.reciVecCoord(:, 2) = reciYcoord(:); 
    S.reciVecCoord(:, 3) = reciZcoord(:); 
    vecLength = vecLength3D(:); % get the length of these reciprocal lattice vectors
%% spline interpolation and FFT to get thetas (p(q)*rho) in reciprocal space
    % D2yDx2 of spline functions are generated in vdWDFinitialize_InputKernel
    [ps, DpDq0s] = spline_interpolation_more(S.vdWDF_qmesh, S.vdWDF_q0, S.vdWDF_D2yDx2); % solve p
    thetas = zeros(nnr, qnum);
    S.vdWDF_thetasFFT = zeros(nnr, qnum);
    for q = 1:qnum
        thetas(:, q) = ps(:, q).*rho(:, 1);
        thetasFFT = fftn(reshape(thetas(:, q), S.Nx,S.Ny,S.Nz));
        S.vdWDF_thetasFFT(:, q) = thetasFFT(:);
    end
    S.vdWDF_thetasFFT = S.vdWDF_thetasFFT*(1/nnr); % Omega/nnr is differential element
%% find the unique reciprocal lattice vector lengthes, then find the value of kernel functions at them
    [uniVecLen, numbers, indexes] = unique(vecLength); % remove repeat lengths
    S.vdWuniReciVecLength = uniVecLen; % to be used in calculating stress
    S.vdWuniReciVecIndexes = indexes; % to be used in calculating stress
    numUniVecLen = size(uniVecLen, 1);
    boolLargerLimit = uniVecLen > S.vdWDF_dk*S.vdWDF_Nrpoints;
    if ismember(1, boolLargerLimit)
        fprintf('in vdWDF_energy, there are reciprocal lattice vectors whose length are larger than limit. Please try using a coarser mesh.\n');
    end
    uniVecLenLLimit = uniVecLen(~boolLargerLimit); % actually, if there is reciprocal vector longer than the length linit,
    % the computation should be stopped and a larger kernel functions
    % should be generated.
    kernelofLengths = interpolate_kernel(uniVecLenLLimit, S.vdWDF_kernel, S.vdWDF_d2Phidk2, S.vdWDF_dk);
    kernelofAllLengths = zeros(numUniVecLen, qnum, qnum);
    kernelofAllLengths(1:size(kernelofLengths, 1), :, :) = kernelofLengths(:, :, :);
    kernelOnPoints = kernelofAllLengths(indexes(:), :, :); % all kernel function values on all reci lattice vectors
    thetaOnPoints = S.vdWDF_thetasFFT;
    u = zeros(ngm, qnum);
    for q2 = 1:qnum
        for q1 = 1:qnum
            u(:, q2) = u(:, q2) + kernelOnPoints(:, q2, q1).*thetaOnPoints(:, q1);
        end
        S.vdWenergy = S.vdWenergy + thetaOnPoints(:, q2)'*u(:, q2); % dot multiplication of complex vector contains conjugate
    end
%% modify scaling factors: (2pi^3/Omega)[differential element]*Omega*Omega[from theta after FFT]*(1/(2pi^3))[kernel FT-IFT, make sure correct scale]
    S.vdWenergy = real(S.vdWenergy) * 0.5 * omega;
    S.vdW_u = u; % prepare for solving potential of vdW. theta_beta.*Phi_alphabetaqnum = size(S.vdWDF_qmesh, 1);
    
    for q = 1:qnum % transfer u back to real space from reciprocal space
        uReci3D = reshape(S.vdW_u(:, q), S.Nx, S.Ny, S.Nz);
        u3D = ifftn(uReci3D)*nnr; % the position of const 1/(Nx*Ny*Nz) is in ifftn in MATLAB 
        %but the const is in fwfft(cfft3d R->G) in FORTRAN
        S.vdW_u(:, q) = u3D(:);
    end
    
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

function kernelofLnegths = interpolate_kernel(uniVecLen, kernel, D2PhiDk2, dk)
    timeofdk = floor(uniVecLen / dk) + 1;
    uniLenNum = size(uniVecLen, 1);
    qnum = size(kernel, 2);
    A = (dk*((timeofdk - 1) + 1.0) - uniVecLen)/dk;
    B = (uniVecLen - dk*(timeofdk - 1))/dk;
    C = (A.^3 - A)*dk^2 / 6.0;
    D = (B.^3 - B)*dk^2 / 6.0;
    kernelofLnegths = zeros(uniLenNum, qnum, qnum);
    for q1 = 1:qnum
        for q2 = 1:q1
            kernelofLnegths(:, q1, q2) = A.*kernel(timeofdk(:), q1, q2) + B.*kernel(timeofdk(:)+1, q1, q2)...
                + (C.*D2PhiDk2(timeofdk(:), q1, q2) + D.*D2PhiDk2(timeofdk(:)+1, q1, q2));
            kernelofLnegths(:, q2, q1) = kernelofLnegths(:, q1, q2);
        end
    end
end