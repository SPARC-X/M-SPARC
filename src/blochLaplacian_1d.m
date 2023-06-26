function [DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,kptvec)
% @ brief    Calculates each component of laplacian in 1D
% @ authors
%         Abhiraj Sharma <asharma424@gatech.edu>
%         Qimen Xu <qimenxu@gatech.edu>
%         Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param kptvec      k-point vector for the current Block diagonalized problem
% @param DLii        Discrete laplacian component in 1D along ith direction
% @param DGi         Discrete gradient component in 1D along ith direction
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%============================================================================
Nx = S.Nx; Ny = S.Ny; Nz = S.Nz;

% Phase factors
if kptvec(1) == 0
	x_phasefac_l = 1.0;
	x_phasefac_r = 1.0;
else
	x_phasefac_l = exp(-1i*kptvec(1)*S.L1);
	x_phasefac_r = exp(1i*kptvec(1)*S.L1);
end

if kptvec(2) == 0
	y_phasefac_l = 1.0;
	y_phasefac_r = 1.0;
else
	y_phasefac_l = exp(-1i*kptvec(2)*S.L2);
	y_phasefac_r = exp(1i*kptvec(2)*S.L2);
end

if kptvec(3) == 0
	z_phasefac_l = 1.0;
	z_phasefac_r = 1.0;
else
	z_phasefac_l = exp(-1i*kptvec(3)*S.L3);
	z_phasefac_r = exp(1i*kptvec(3)*S.L3);
end


% D_xx laplacian in 1D
%-----------------------
V = S.V_11;

if S.BCx == 0
	V(S.isOutl_11) = V(S.isOutl_11) * x_phasefac_l;
	V(S.isOutr_11) = V(S.isOutr_11) * x_phasefac_r;
end

% Create discretized Laplacian
DL11 = sparse(S.I_11,S.II_11,V,Nx,Nx);

% D_yy laplacian in 1D
%-----------------------
V = S.V_22;

if S.BCy == 0
	V(S.isOutl_22) = V(S.isOutl_22) * y_phasefac_l;
	V(S.isOutr_22) = V(S.isOutr_22) * y_phasefac_r;  
end

% Create discretized Laplacian
DL22 = sparse(S.I_22,S.II_22,V,Ny,Ny);

% D_zz laplacian in 1D
%-----------------------

V = S.V_33;

if S.BCz == 0
	V(S.isOutl_33) = V(S.isOutl_33) * z_phasefac_l;
	V(S.isOutr_33) = V(S.isOutr_33) * z_phasefac_r;  
end

% Create discretized Laplacian
DL33 = sparse(S.I_33,S.II_33,V,Nz,Nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DG1 = zeros(Nx,Nx);
DG2 = zeros(Ny,Ny);
DG3 = zeros(Nz,Nz);

if S.cell_typ > 1
	% x-direction
	%-------------

	V = S.V_1;

	if S.BCx == 0
		V(S.isOutl_1) = V(S.isOutl_1) * x_phasefac_l;
		V(S.isOutr_1) = V(S.isOutr_1) * x_phasefac_r;  
	end

	% Create discretized Laplacian
	DG1 = sparse(S.I_1,S.II_1,V,Nx,Nx);

	% y-direction
	%-------------

	V = S.V_2;

	if S.BCy == 0
		V(S.isOutl_2) = V(S.isOutl_2) * y_phasefac_l;
		V(S.isOutr_2) = V(S.isOutr_2) * y_phasefac_r;  
	end

	% Create discretized Laplacian
	DG2 = sparse(S.I_2,S.II_2,V,Ny,Ny);

	% z-direction
	%-------------

	V = S.V_3;

	if S.BCz == 0
		V(S.isOutl_3) = V(S.isOutl_3) * z_phasefac_l;
		V(S.isOutr_3) = V(S.isOutr_3) * z_phasefac_r;  
	end

	% Create discretized Laplacian
	DG3 = sparse(S.I_3,S.II_3,V,Nz,Nz);
	
end
