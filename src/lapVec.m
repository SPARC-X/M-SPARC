function DLX = lapVec(DL11,DL22,DL33,DG1,DG2,DG3,X,S)
% @ brief    Calculates laplacian vector product using 
%            Kronecker product method
% @ authors
%         Abhiraj Sharma <asharma424@gatech.edu>
%         Qimen Xu <qimenxu@gatech.edu>
%         Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param DLii        Discrete laplacian component in 1D along ith direction
% @param DGi         Discrete gradient component in 1D along ith direction
% @param X           N x Ns matrix of eigenfunctions of Hamiltonian
% @param DLX         Discrete laplacian times X
%
% @ references
%              "On real-space Density Functional Theory for 
%               non-orthogonal crystal systems: Kronecker product 
%               formulation of the kinetic energy operator (Sharma et. al. 2018)"
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%=====================================================================================
Nx = S.Nx;Ny = S.Ny;Nz = S.Nz;
if S.cell_typ == 1 % orthogonal systems
	X1 = reshape(X,Nx,Ny,[]);
	l_zs = size(X1,3);
	Hlx1 = zeros(Nx,Ny,l_zs);
	for i=1:l_zs
		Hlx1(:,:,i) = DL11*X1(:,:,i) + X1(:,:,i)*DL22.';
	end
	Hlx1 = reshape(Hlx1,Nx*Ny,Nz,[]) ;
	X2 = reshape(X,Nx*Ny,Nz,[]);
	l_s = size(X2,3);
	Hlx2 = zeros(Nx*Ny,Nz,l_s);
	for i=1:l_s
		Hlx2(:,:,i) = X2(:,:,i)*DL33.';
	end    
	DLX = reshape((Hlx1 + Hlx2),[],l_s);
elseif S.cell_typ == 2
	T11 = S.lapc_T(1,1);T22 = S.lapc_T(2,2);T33 = S.lapc_T(3,3);T12 = S.lapc_T(1,2);
	T23 = S.lapc_T(2,3);T13 = S.lapc_T(1,3);
	
	X1 = reshape(X,Nx,Ny,[]);
	l_zs = size(X1,3);
	Hlx1 = zeros(Nx,Ny,l_zs);
	Hlx2 = zeros(Nx,Ny,l_zs);
	for i=1:l_zs
		Hlx1(:,:,i) = T11*DL11*X1(:,:,i) + T22*X1(:,:,i)*DL22.' + T12*DG1*X1(:,:,i)*DG2.' ;
		Hlx2(:,:,i) = T13*DG1*X1(:,:,i) + T23*X1(:,:,i)*DG2.' ;
	end
	Hlx1 = reshape(Hlx1,Nx*Ny,Nz,[]);
	Hlx2 = reshape(Hlx2,Nx*Ny,Nz,[]);

	X2 = reshape(X,Nx*Ny,Nz,[]);
	l_s = size(X2,3);
	for i=1:l_s
		Hlx2(:,:,i) = T33*X2(:,:,i)*DL33.' + Hlx2(:,:,i)*DG3.';
	end
	
	DLX = reshape((Hlx1 + Hlx2),[],l_s);    
elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
	DLX = lapVec_cychel(DL11,DL22,DL33,DG1,DG2,DG3,X,S);
end
	
	

