function [S] = vdWDF_Initial_GenKernel(S)
% @file    vdWDF_Initial_GenKernel.m
% @brief   This file contains the functions for generating the needed model 
%          kernel functions in reciprocal space and the value of spline
%          functions at model energy ratios.
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
%% Initialization, set parameters and grids
    S.vdWDF_Nrpoints = 1024; %% radial points for composing Phi functions in real and reciprocal space
    S.vdWDF_rMax = 100.0; %% max radius in real space and minimum k point 2*pi/r_max in reciprocal space
    S.vdWDF_dr = S.vdWDF_rMax/S.vdWDF_Nrpoints;
    S.vdWDF_dk = 2.0*pi/S.vdWDF_rMax;
    S.vdWDF_Nqs = 20;
    S.vdWDF_qmesh = [
    1.0e-5            0.0449420825586261 0.0975593700991365 0.159162633466142
    0.231286496836006 0.315727667369529  0.414589693721418  0.530335368404141
    0.665848079422965 0.824503639537924  1.010254382520950  1.227727621364570
    1.482340921174910 1.780437058359530  2.129442028133640  2.538050036534580
    3.016440085356680 3.576529545442460  4.232271035198720  5.0];
    S.vdWDF_qmesh = reshape(S.vdWDF_qmesh', [], 1);
    S.vdWDF_kernel = zeros(1 + S.vdWDF_Nrpoints, S.vdWDF_Nqs, S.vdWDF_Nqs); %% kernal Phi, index 0, reciprocal
    S.vdWDF_d2Phidk2 = zeros(1 + S.vdWDF_Nrpoints, S.vdWDF_Nqs, S.vdWDF_Nqs); %% 2nd derivative of kernal
    S.vdWenergy = 0.0;
    fprintf('\n Begin generating kernel functions of vdW-DF...\n');
%% Generate Wab(a, b) function
    NintegratePoints = 256; % Number of imagine frequency points integrated for building kernel
    aMin = 0.0; aMax = 64.0; % scope of integration of imagine frequency (a, b)
    aPoints = zeros(NintegratePoints, 1); % imagine frequency points, column
    [weights, aPoints] = prepGaussQuad(NintegratePoints, aMin, aMax, aPoints); % Gaussian quadrature integration points and weights
    aPoints = tan(aPoints);
    aPoints2 = aPoints(:).*aPoints(:);
    weights(:) = weights(:).*(1 + aPoints2(:));
    cos_a = cos(aPoints);
    sin_a = sin(aPoints);
    Wab = 2.0 * (weights*weights') ./ (aPoints*aPoints'); %% Wab(NintegratePoints, NintegratePoints)
    part1 = ((3.0 - aPoints2)*aPoints') .* (sin_a*cos_a');
    part2 = (aPoints*(3.0 - aPoints2)') .* (cos_a*sin_a');
    part3 = (repmat(aPoints2, 1, NintegratePoints) + repmat(aPoints2', NintegratePoints, 1) - 3.0) .* (sin_a*sin_a');
    part4 = -3.0 * (aPoints*aPoints') .* (cos_a*cos_a');
    Wab = Wab .* (part1 + part2 + part3 + part4); %% attention: Wab is diff from W_ab calculated by QE!
    % W_ab from QE is not accurate
%% Compute kernal function Phi in reciprocal space
    for alphaQ = 1:S.vdWDF_Nqs
        for betaQ = 1:alphaQ
            fprintf('Generating kernel function K(p %d, p %d)\n', alphaQ, betaQ);
            d1 = (1:S.vdWDF_Nrpoints)' * S.vdWDF_dr * S.vdWDF_qmesh(alphaQ);
            d2 = (1:S.vdWDF_Nrpoints)' * S.vdWDF_dr * S.vdWDF_qmesh(betaQ);
            phi = phiValue(d1, d2, Wab, aPoints, aPoints2);
            reciPhi = radialFFT(phi, S.vdWDF_dr, S.vdWDF_dk);
            S.vdWDF_kernel(:, alphaQ, betaQ) = reciPhi(:);
            S.vdWDF_kernel(:, betaQ, alphaQ) = reciPhi(:);
            d2reciPhi = d2ForSplines(reciPhi, S.vdWDF_dk); %% set_up_splines
            S.vdWDF_d2Phidk2(:, alphaQ, betaQ) = d2reciPhi(:);
            S.vdWDF_d2Phidk2(:, betaQ, alphaQ) = d2reciPhi(:);
        end
    end
%% generate D2yDx2 for future spline interpolation in getQ0onGrid
    S.vdWDF_D2yDx2 = splineFunc_2Deri_atQmesh(S.vdWDF_qmesh); % initialize_spline_interpolation in QE

end

function [weights, aPoints] = prepGaussQuad(NintegratePoints, aMin, aMax, aPoints)
% should be a column
    weights = zeros(NintegratePoints, 1);
    Npoints = (NintegratePoints + 1)/2;
    midPoint = 0.5 * (atan(aMin) + atan(aMax));
    lengthScope = 0.5 * (atan(aMax) - atan(aMin));
    rootVec = (1:Npoints)';
    rootVec = cos((pi*(rootVec - 0.25) / (NintegratePoints + 0.5)));
    for iPoint = 1:Npoints
        root = rootVec(iPoint);
        rootFlag = 1;
        while rootFlag %% vectorized, for all Npoints
            poly1 = 1.0;
            poly2 = 0.0;
            for iPoly = 1:NintegratePoints
                poly3 = poly2;
                poly2 = poly1;
                poly1 = ((2.0*iPoly - 1.0)*root*poly2 - (iPoly - 1.0)*poly3) / iPoly;
            end
            dpdx = NintegratePoints * (root*poly1 - poly2) / (root^2 - 1.0);
            last_root = root;
            root = last_root - poly1/dpdx;
            if (abs(root - last_root) <= 1e-14)
                rootFlag = 0;
            end
        end
        aPoints(iPoint) = midPoint - lengthScope*root;
        aPoints(NintegratePoints + 1 - iPoint) = midPoint + lengthScope*root;
        weights(iPoint) = 2.0*lengthScope / ((1.0-root^2) * dpdx^2);
        weights(NintegratePoints + 1 - iPoint) = weights(iPoint);
    end
end

function phi = phiValue(d1, d2, Wab, aPoints, aPoints2)
    Nrpoints = size(d1, 1);
    NintegratePoints = size(aPoints, 1);
    phi = zeros(Nrpoints, 1);
    for r = 1:Nrpoints
        nu = aPoints2 ./ (2.0 * hFunction(aPoints ./ d1(r)));
        nu1 = aPoints2 ./ (2.0 * hFunction(aPoints ./ d2(r)));
        % wVec = nu(a_1); yVec = nu1(a_i); xVec = nu(b_i); zVec = nu1(b_i);
        % a_i: row index; b_i: column index
        wPxMatrix = repmat(nu, 1, NintegratePoints) + repmat(nu', NintegratePoints, 1);
        yPzMatrix = repmat(nu1, 1, NintegratePoints) + repmat(nu1', NintegratePoints, 1);
        wPyMatrix = repmat((nu + nu1), 1, NintegratePoints);
        xPzMatrix = repmat((nu + nu1)', NintegratePoints, 1);
        wPzMatrix = repmat(nu, 1, NintegratePoints) + repmat(nu1', NintegratePoints, 1);
        yPxMatrix = repmat(nu1, 1, NintegratePoints) + repmat(nu', NintegratePoints, 1); 
        Tmatrix = (1.0./wPxMatrix + 1.0./yPzMatrix) .* (1.0./(wPyMatrix.*xPzMatrix) + 1.0./(wPzMatrix.*yPxMatrix));
        phi(r) = sum(sum(Tmatrix .* Wab));
    end
    phi = 1.0 / pi^2 * phi;
end

function reciPhi = radialFFT(phi, dr, dk)
    Nrpoints = size(phi, 1);
    reciPhi = zeros(Nrpoints + 1, 1); %% index from 0 to Nrpoints
    r = (1:Nrpoints)' * dr;
    reciPhi(1) = 4*pi * sum(phi .* (r.^2));
    reciPhi(1) = reciPhi(1) - 4*pi *0.5*r(Nrpoints)^2*phi(Nrpoints);
    k = (1:Nrpoints)' * dk;
    for indexK = 2:Nrpoints + 1
        besselKmulR = 4*pi * sin(k(indexK - 1)*r) / k(indexK - 1); %% spherical bessel function * radius
        reciPhi(indexK) = sum(phi .* r .* besselKmulR);
        reciPhi(indexK) = reciPhi(indexK) - 0.5*phi(Nrpoints)*r(Nrpoints)*besselKmulR(Nrpoints);
    end
    reciPhi = reciPhi * dr;
end

function d2reciPhi = d2ForSplines(reciPhi, dk)
    Nrpoints = size(reciPhi, 1) - 1;
    d2reciPhi = zeros(Nrpoints + 1, 1);
    tempArray = zeros(Nrpoints + 1, 1);
    for indexK = 2:Nrpoints
        temp1 = 0.5;
        temp2 = temp1*d2reciPhi(indexK - 1) + 2.0;
        d2reciPhi(indexK) = (temp1 - 1.0) / temp2;
        tempArray(indexK) = (reciPhi(indexK + 1) - reciPhi(indexK))/dk - (reciPhi(indexK) - reciPhi(indexK - 1))/dk;
        tempArray(indexK) = (6.0*tempArray(indexK)/(2*dk) - temp1*tempArray(indexK - 1)) / temp2;
    end
    d2reciPhi(Nrpoints + 1) = 0.0;
    for indexK = Nrpoints:-1:1
        d2reciPhi(indexK) = d2reciPhi(indexK)*d2reciPhi(indexK + 1) + tempArray(indexK);
    end
end

function hResultArray = hFunction(aPointsDivideD)
    g1 = 4.0*pi/9.0;
    hResultArray = 1.0 - exp(-g1*aPointsDivideD.^2);
end


function D2yDx2 = splineFunc_2Deri_atQmesh(qmesh)
    qnum = size(qmesh, 1);
    D2yDx2 = zeros(qnum);
    for q = 1:qnum
        y = zeros(qnum, 1);
        y(q) = 1;
        tempArray = zeros(qnum, 1);
        for q2 = 2:qnum - 1
            temp1 = (qmesh(q2) - qmesh(q2-1))/(qmesh(q2+1) - qmesh(q2-1));
            temp2 = temp1*D2yDx2(q, q2 - 1) + 2.0;
            D2yDx2(q, q2) = (temp1 - 1.0)/temp2;
            tempArray(q2) = (y(q2+1) - y(q2))/(qmesh(q2+1) - qmesh(q2))...
             - (y(q2) - y(q2-1))/(qmesh(q2) - qmesh(q2-1));
            tempArray(q2) = (6.0*tempArray(q2)/(qmesh(q2+1) - qmesh(q2-1))...
             - temp1*tempArray(q2-1))/temp2;
        end
        for q2 = qnum - 1:-1:1
            D2yDx2(q, q2) = D2yDx2(q, q2)*D2yDx2(q, q2+1) + tempArray(q2);
        end
    end
end