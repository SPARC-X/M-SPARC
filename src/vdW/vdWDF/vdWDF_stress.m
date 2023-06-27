function S = vdWDF_stress(S)
% @file    vdWDF_stress.m
% @brief   This file contains the functions computing vdWDF stress term
% @authors Boqin Zhang <bzhang376@gatech.edu>
%          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% Reference:
% Dion, Max, Henrik Rydberg, Elsebeth Schröder, David C. Langreth, and Bengt I. Lundqvist. 
% "Van der Waals density functional for general geometries." 
% Physical review letters 92, no. 24 (2004): 246401.
% Román-Pérez, Guillermo, and José M. Soler. 
% "Efficient implementation of a van der Waals density functional: application to double-wall carbon nanotubes." 
% Physical review letters 103, no. 9 (2009): 096102.
% Lee, Kyuho, Éamonn D. Murray, Lingzhu Kong, Bengt I. Lundqvist, and David C. Langreth. 
% "Higher-accuracy van der Waals density functional." Physical Review B 82, no. 8 (2010): 081101.
% Thonhauser, T., S. Zuluaga, C. A. Arter, K. Berland, E. Schröder, and P. Hyldgaard. 
% "Spin signature of nonlocal correlation binding in metal-organic frameworks." 
% Physical review letters 115, no. 13 (2015): 136402.
% Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
% ==============================================================================================
    S.vdWstress = zeros(3);
    if S.nspin == 1
        vdWstrGradTerm = vdWDF_stress_gradient(S);
    else
        vdWstrGradTerm = SvdWDF_stress_gradient(S);
    end
    vdWstrKernelTerm = vdWDF_stress_kernel(S);
    S.vdWstress = vdWstrGradTerm + vdWstrKernelTerm; % initial sign is minus, which is wrong
    for l=1:3
        for m = 1:l-1
            S.vdWstress(m, l) = S.vdWstress(l, m);
        end
    end
    S.vdWstress = real(S.vdWstress);
end

function stressGrad = vdWDF_stress_gradient(S)
    stressGrad = zeros(3);
    DpDq0s = S.vdW_DpDq0s;
    nnr = size(DpDq0s, 1);
    qnum = size(DpDq0s, 2);
    gradRhoLength = sqrt(S.Drho(:, 1).^2 + S.Drho(:, 2).^2 + S.Drho(:, 3).^2);
    u = S.vdW_u;
    Dq0Dgradrho = S.vdWDF_Dq0Dgradrho;
    prefactor = zeros(nnr, 1);
    for q = 1:qnum
        prefactor = prefactor + u(:, q).*DpDq0s(:, q).*Dq0Dgradrho(:);
    end
    gradLarger0 = gradRhoLength > 0.0;
    prefactorGradLarger0 = prefactor(gradLarger0)./gradRhoLength(gradLarger0);
    for l=1:3
        for m=1:l
            productlm = S.Drho(gradLarger0, l).*S.Drho(gradLarger0, m);
            stressGrad(l, m) = stressGrad(l, m) - sum(prefactorGradLarger0.*productlm);
        end
    end
    stressGrad = stressGrad/nnr;
end

function stressGrad = SvdWDF_stress_gradient(S)
    stressGrad = zeros(3);
    DpDq0s = S.vdW_DpDq0s;
    nnr = size(DpDq0s, 1);
    qnum = size(DpDq0s, 2);
    gradRhoLengthUp = vecnorm(S.DrhoUp, 2, 2);
    gradRhoLengthDn = vecnorm(S.DrhoDn, 2, 2);
    gradLarger0 = (gradRhoLengthUp > 0.0) & (gradRhoLengthDn > 0.0);
    u = S.vdW_u;
    Dq0Dgradrho = S.vdWDF_Dq0Dgradrho;
    prefactor = zeros(nnr, 2);
    for q = 1:qnum
        prefactor = prefactor + (u(:, q).*DpDq0s(:, q)) .* Dq0Dgradrho; % col 1: prefactor_up; col 2: prefactor_dn
    end
    prefactorGradLarger0 = prefactor(gradLarger0, :)./[gradRhoLengthUp(gradLarger0), gradRhoLengthDn(gradLarger0)];
    for l=1:3
        for m=1:l
            productlm = [S.DrhoUp(gradLarger0, l).*S.DrhoUp(gradLarger0, m), S.DrhoDn(gradLarger0, l).*S.DrhoDn(gradLarger0, m)];
            stressGrad(l, m) = stressGrad(l, m) - sum(sum(prefactorGradLarger0.*productlm));
        end
    end
    stressGrad = stressGrad/nnr;
end


function stressKernel = vdWDF_stress_kernel(S)
    qnum = size(S.vdW_DpDq0s, 2);
    uniReciVecLen = S.vdWuniReciVecLength;
    uniReciVecIndexes = S.vdWuniReciVecIndexes;
    numUniVecLen = size(uniReciVecLen, 1);
    boolLargerLimit = uniReciVecLen > S.vdWDF_dk*S.vdWDF_Nrpoints; % find reciprocal vectors length less than limit of spline interpolation
    if ismember(1, boolLargerLimit)
        fprintf('in vdWDF_energy, there are reciprocal lattice vectors whose length are larger than limit.\n');
    end
    uniVecLenLLimit = uniReciVecLen(~boolLargerLimit);
    DkernelDkofLengths = interpolate_DkernelDk(uniVecLenLLimit, S.vdWDF_kernel, S.vdWDF_d2Phidk2, S.vdWDF_dk);
    DkernelDkofAllLengths = zeros(numUniVecLen, qnum, qnum);
    DkernelDkofAllLengths(1:size(DkernelDkofLengths, 1), :, :) = DkernelDkofLengths(:, :, :);
    DkernelDkOnPoints = DkernelDkofAllLengths(uniReciVecIndexes(:), :, :); % all kernel function values on all reci lattice vectors
    thetaOnPoints = S.vdWDF_thetasFFT;
    reciVecCoord = S.reciVecCoord; % coordinate of g vectors
    ngm = size(reciVecCoord, 1);
    reciVecLen = uniReciVecLen(uniReciVecIndexes(:)); % length = 0 needs to be excluded
    reciVecLenLZero = reciVecLen > 1e-8;
    
    theta_DkernelDk_conjTheta = zeros(ngm, 1);
    for q2 = 1:qnum
        for q1 = 1:qnum
            theta_DkernelDk_conjTheta = theta_DkernelDk_conjTheta + ...
                thetaOnPoints(:, q1).*DkernelDkOnPoints(:,q1,q2).*conj(thetaOnPoints(:, q2));
        end
    end
    stressKernel = zeros(3);
    for l = 1:3
        for m = 1:l
            stressKernel(l, m) = -0.5*sum(theta_DkernelDk_conjTheta(reciVecLenLZero)...
                .*reciVecCoord(reciVecLenLZero, l).*reciVecCoord(reciVecLenLZero, m)./reciVecLen(reciVecLenLZero));
        end
    end
end

function DkernelDkofLengths = interpolate_DkernelDk(uniVecLen, kernel, D2PhiDk2, dk)
    timeofdk = floor(uniVecLen / dk) + 1;
    uniLenNum = size(uniVecLen, 1);
    qnum = size(kernel, 2);
    A = (dk*((timeofdk - 1) + 1.0) - uniVecLen)/dk;
    B = (uniVecLen - dk*(timeofdk - 1))/dk;
    dAdk = -1.0/dk;
    dBdk = 1.0/dk;
    dCdk = -((3*A.^2 - 1.0)/6.0)*dk;
    dDdk = ((3*B.^2 -1.0)/6.0)*dk;
    DkernelDkofLengths = zeros(uniLenNum, qnum, qnum);
    for q1 = 1:qnum
        for q2 = 1:q1
            DkernelDkofLengths(:, q1, q2) = dAdk*kernel(timeofdk(:), q1, q2) + dBdk*kernel(timeofdk(:)+1, q1, q2)...
                + (dCdk.*D2PhiDk2(timeofdk(:), q1, q2) + dDdk.*D2PhiDk2(timeofdk(:)+1, q1, q2));
            DkernelDkofLengths(:, q2, q1) = DkernelDkofLengths(:, q1, q2);
        end
    end
end