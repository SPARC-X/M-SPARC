function S = initialization(filename)
% @brief    initialization(filename) reads data from input file 'filename'
%           and creates a struct 'S' to store all the data read and initialized.
%
% @param filename   The input filename.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech

% Set up inpt defaults
S = inpt_defaults();

% Read .inpt file
S = read_inpt(S, filename);

% Read .ion file
S = read_ion(S, filename);

% Read pseudopotential files
for ityp = 1:S.n_typ
	if (S.Atm(ityp).psptyp == 0)
		% WARNING: for TM format, lloc has to be set before reading
		S = ReadPseudoPot(S.Atm(ityp).lloc,S.Atm(ityp).typ);
	elseif (S.Atm(ityp).psptyp == 1)
		S = readPseudopot(S, ityp, S.Atm(ityp).psdfname, S.Atm(ityp).typ);
	else
		error('Cannot recognize pseudopotential format!');
	end
end

% Set up more defaults based on input files
S = setup_defaults(S, filename);

% Calculate rb
S = Calculate_rb(S);

% write initialized parameters into output file
S = Write_output_init(S, filename);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = Calculate_rb(S)
% Starting and ending indices of b-region
if (S.cell_typ == 1 || S.cell_typ == 2)
	pos_atm_x = 0; % atom location in x-direction
	pos_atm_y = 0; % atom location in y-direction
	pos_atm_z = 0; % atom location in z-direction
	rb_up_x = 12;
	f_rby = @(y) y;
	rb_up_y = f_rby(12);
	rb_up_z = 12;
	
elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
	pos_atm_x = S.xmax_at; % maximum R coordinate of any atom
	pos_atm_y = 0; % atom location in theta-direction
	pos_atm_z = 0; % atom location in z-direction
	rb_up_x = S.xvac; % Radial direction vacuum
	f_rby = @(y) acos(1 - y^2/(2*pos_atm_x^2));
	rb_up_y = f_rby(12); % Theta direction
	rb_up_z = 12; % z-direction
	
end
ii_s_temp = -ceil(rb_up_x/S.dx);
ii_e_temp = ceil(rb_up_x/S.dx);
jj_s_temp = -ceil(rb_up_y/S.dy);
jj_e_temp = ceil(rb_up_y/S.dy);
kk_s_temp = 0;
kk_e_temp = ceil(rb_up_z/S.dz);
xx_temp = pos_atm_x + (ii_s_temp-S.FDn:ii_e_temp+S.FDn)*S.dx;
yy_temp = pos_atm_y + (jj_s_temp-S.FDn:jj_e_temp+S.FDn)*S.dy;
zz_temp = pos_atm_z + (kk_s_temp-S.FDn:kk_e_temp+S.FDn)*S.dz;
[XX_3D_temp,YY_3D_temp,ZZ_3D_temp] = ndgrid(xx_temp,yy_temp,zz_temp);
Nx = (ii_e_temp-ii_s_temp)+1;
Ny = (jj_e_temp-jj_s_temp)+1;
Nz = (kk_e_temp-kk_s_temp)+1;
% Find distances
dd_temp = calculateDistance(XX_3D_temp,YY_3D_temp,ZZ_3D_temp,pos_atm_x,pos_atm_y,pos_atm_z,S);

% Find integration weights
W_temp = IntgWts(Nx,Ny,Nz,1,1,1,xx_temp(S.FDn+1),S); % 1 - dirichlet BC on the boundary nodes
W_temp = reshape(W_temp,Nx,Ny,Nz);

% Find VJ and bJ
for ityp = 1:S.n_typ
	V_PS_temp = zeros(size(dd_temp));
	IsLargeThanRmax = dd_temp > S.Atm(ityp).r_grid_vloc(end);
	V_PS_temp(IsLargeThanRmax) = -S.Atm(ityp).Z;
	V_PS_temp(~IsLargeThanRmax) = interp1(S.Atm(ityp).r_grid_vloc, ...
		S.Atm(ityp).r_grid_vloc.*S.Atm(ityp).Vloc, dd_temp(~IsLargeThanRmax), 'spline');
	
	V_PS_temp = V_PS_temp./dd_temp;
	V_PS_temp(dd_temp<S.Atm(ityp).r_grid_vloc(2)) = S.Atm(ityp).Vloc(1);
	II_temp = 1+S.FDn : size(V_PS_temp,1)-S.FDn;
	JJ_temp = 1+S.FDn : size(V_PS_temp,2)-S.FDn;
	KK_temp = 1+S.FDn : size(V_PS_temp,3)-S.FDn;
	
	b_temp = pseudochargeDensity_atom(V_PS_temp,II_temp,JJ_temp,KK_temp,xx_temp(1),S);
	b_temp = -b_temp / (4*pi);
	err_rb = 100;
	count = 1;
	rb_x = S.Atm(ityp).rc;
	rb_y = f_rby(S.Atm(ityp).rc);
	rb_z = S.Atm(ityp).rc;
    rb_x = ceil(rb_x/S.dx-1e-12)*S.dx;
    rb_y = ceil(rb_y/S.dy-1e-12)*S.dy;
    rb_z = ceil(rb_z/S.dz-1e-12)*S.dz;
	fprintf(' Finding rb for %s ...\n',S.Atm(ityp).typ);
	while (err_rb > S.pseudocharge_tol && count <= 100 && rb_x <= rb_up_x && rb_y <= rb_up_y && rb_z <= rb_up_z )
		rb_x = rb_x + S.dx;
		rb_z = rb_z + S.dz;
		rb_y = f_rby(max(rb_x,rb_z));
		ii_rb = -1*ii_s_temp+S.FDn-floor(rb_x/S.dx)+1:-1*ii_s_temp+S.FDn+floor(rb_x/S.dx)+1;
		jj_rb = -1*jj_s_temp+S.FDn-floor(rb_y/S.dy)+1:-1*jj_s_temp+S.FDn+floor(rb_y/S.dy)+1;
		kk_rb = S.FDn+1:S.FDn+floor(rb_z/S.dz)+1; 
		err_rb = abs(sum(sum(sum(W_temp(ii_rb-S.FDn,jj_rb-S.FDn,kk_rb-S.FDn).*b_temp(ii_rb,jj_rb,kk_rb))))*2 + S.Atm(ityp).Z);
		fprintf(' rb = {%.3f %.3f %.3f}, int_b = %.15f, err_rb = %.3e\n',rb_x,rb_y,rb_z,2*sum(sum(sum(W_temp(ii_rb-S.FDn,jj_rb-S.FDn,kk_rb-S.FDn).*b_temp(ii_rb,jj_rb,kk_rb)))),err_rb);
		count = count + 1;
	end
	
	assert(rb_x<=rb_up_x && rb_y<=rb_up_y && rb_z<=rb_up_z,'Need to increase upper bound for rb!');
	S.Atm(ityp).rb_x = rb_x;
	S.Atm(ityp).rb_y = rb_y;
	S.Atm(ityp).rb_z = rb_z;
	% S.Atm(ityp).rb_x = ceil(rb_x/S.dx-1e-12)*S.dx; % + S.dx;
	% S.Atm(ityp).rb_y = ceil(rb_y/S.dy-1e-12)*S.dy; % + S.dy;
	% S.Atm(ityp).rb_z = ceil(rb_z/S.dz-1e-12)*S.dz; % + S.dz;
    fprintf(' rb = {%.3f %.3f %.3f}\n',S.Atm(ityp).rb_x,S.Atm(ityp).rb_y,S.Atm(ityp).rb_z);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = setup_defaults(S, filename)
% After reading the input files, set up remaining default values
S.temp_tol = 1e-12;

% Density tolerance for exchange-correlation
S.xc_rhotol = 1e-14;

% exchange correlation
if strcmp(S.XC, 'LDA_PW') || strcmp(S.XC, 'LDA')
	S.xc = 0;
elseif strcmp(S.XC, 'LDA_PZ')
	S.xc = 1; % TODO: Implement PZ LDA!
elseif strcmp(S.XC, 'GGA_PBE') || strcmp(S.XC, 'GGA')
	S.xc = 2;
end

% calculate Nelectron
S.Nelectron = 0;
for ityp = 1:S.n_typ
	S.Nelectron = S.Nelectron + S.Atm(ityp).Z * S.Atm(ityp).n_atm_typ;
end

% Save the initial positions, since later S.Atoms will be changed
S.Atoms_init = S.Atoms;

% Cychel parameters
if (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
	% Add the folder containing all the Cychel files to the path
	addpath('../Cychel');
	count_prev = 0;
	count = 0;
	for ityp = 1:S.n_typ
		if(S.IsFrac(ityp) == 0)
			S.Atm(ityp).coords = coordinateTransformation_cychel(S,S.Atm(ityp).coords,'cart2noncart_dis');
			count = count + S.Atm(ityp).n_atm_typ;
			S.Atoms(count_prev+1:count,:) = S.Atm(ityp).coords;
			count_prev = count;
		else
			count = count + S.Atm(ityp).n_atm_typ;
			count_prev = count;
		end
	end

	% Find minimum and maximum radial coordinate among all atoms
	S.xmin_at = min(S.Atoms(:,1));
	S.xmax_at = max(S.Atoms(:,1));

	% Find vacuum
	S.xvac = (S.L1 - (S.xmax_at - S.xmin_at))/2;

	% Find inner and outer radius
	S.xin  = S.xmin_at - S.xvac;
	S.xout = S.xmax_at + S.xvac;
	fprintf('Cychel radial direction parameters: \n');
	fprintf('Vacuum %f, Inner radius %f, Outer radius %f\n',S.xvac,S.xin,S.xout);

	% Rotational matrix
	theta1 = S.L2;
	S.RotM1 = [cos(theta1),-sin(theta1),0; sin(theta1),cos(theta1),0; 0 0 1];
	theta2 = S.alph*S.L3;
	S.RotM2 = [cos(theta2),-sin(theta2),0; sin(theta2),cos(theta2),0; 0 0 1];
end

% Check the cell typ (orthogonal or non-orthogonal)
if(abs(dot(S.lat_vec(1,:),S.lat_vec(2,:))) > S.temp_tol ||...
   abs(dot(S.lat_vec(2,:),S.lat_vec(3,:))) > S.temp_tol ||...
   abs(dot(S.lat_vec(3,:),S.lat_vec(1,:))) > S.temp_tol)
	S.cell_typ = 2;
end

S.lat_uvec(1,:) = S.lat_vec(1,:)/norm(S.lat_vec(1,:));
S.lat_uvec(2,:) = S.lat_vec(2,:)/norm(S.lat_vec(2,:));
S.lat_uvec(3,:) = S.lat_vec(3,:)/norm(S.lat_vec(3,:));

% Set up transformation matrices for non orthogonal cells
if S.cell_typ == 2
	% Jacobian
	S.Jacb = det(S.lat_uvec');
	assert(S.Jacb > 0.0,'Volume is negative!');

	% metric_T, Gradient and laplacian transformation matrices
	S.metric_T = S.lat_uvec * S.lat_uvec' ;
	S.metric_T(1,2) = 2*S.metric_T(1,2); 
	S.metric_T(2,3) = 2*S.metric_T(2,3); 
	S.metric_T(1,3) = 2*S.metric_T(1,3);
	S.grad_T = inv(S.lat_uvec') ;
	S.lapc_T = S.grad_T * S.grad_T' ;
	S.lapc_T(1,2) = 2*S.lapc_T(1,2); 
	S.lapc_T(2,3) = 2*S.lapc_T(2,3);
	S.lapc_T(1,3) = 2*S.lapc_T(1,3);

	count_prev = 0;
	count = 0;
	for ityp = 1:S.n_typ
		if(S.IsFrac(ityp) == 0)
			S.Atm(ityp).coords = transpose(S.grad_T * transpose(S.Atm(ityp).coords));
			count = count + S.Atm(ityp).n_atm_typ;
			S.Atoms(count_prev+1:count,:) = S.Atm(ityp).coords;
			count_prev = count;
		else
			count = count + S.Atm(ityp).n_atm_typ;
			count_prev = count;
		end
	end
end

% Set up number of spin
if S.spin_typ == 0
	S.nspin = 1;
else
	S.nspin = 2;
end

% Provide default spin if not provided
if S.nspin == 2
	rng('default');
	for ityp = 1:S.n_typ
		if(S.IsSpin(ityp) == 0)
			S.Atm(ityp).mag = -S.Atm(ityp).Z + 2 * S.Atm(ityp).Z * rand(S.Atm(ityp).n_atm_typ,1);
		end
	end
end

% set up default smearing if not provided
if S.bet < 0
	if S.elec_T_type == 0 % fermi-dirac 
		% The electronic temperature corresponding to 0.1 eV is 1160.452211 K
		S.bet = 27.21138602 / 0.1; % smearing = 0.1 eV = 0.00367493225 Ha, Beta := 1 / smearing
	elseif S.elec_T_type == 1 % gaussian smearing
		% The electronic temperature corresponding to 0.2 eV is 2320.904422 K
		S.bet = 27.21138602 / 0.2; % smearing = 0.2 eV = 0.00734986450 Ha, Beta := 1 / smearing
	end
	S.Temp = 1./(3.166810501187400e-06 * S.bet); 
end

% BCx = 0 -> periodic, BCx = 1 -> dirichlet
if S.BC >= 0    % if user provides BOUNDARY_CONDITION: 1-4
	if(S.BC == 1)
		S.BCx = 1; S.BCy = 1; S.BCz = 1;
	elseif(S.BC == 2)
		S.BCx = 0; S.BCy = 0; S.BCz = 0;
	elseif(S.BC == 3)
		S.BCx = 0; S.BCy = 0; S.BCz = 1;
	elseif(S.BC == 4)
		%S.BCx = 0; S.BCy = 1; S.BCz = 1;
		S.BCx = 1; S.BCy = 1; S.BCz = 0;
	else
		error('Boundary condition should be one among {1,2,3,4}');
	end
elseif S.BCx >= 0 && S.BCy >= 0 && S.BCz >= 0 % if user provides BCx,BCy,BCz
	n_Dirichlet = S.BCx + S.BCy + S.BCz;
	if n_Dirichlet == 0
		S.BC = 2; % Periodic BC in all 3D
	elseif n_Dirichlet == 1
		S.BC = 3; % Surface, Periodic in 2D, Dirichlet in 1D
	elseif n_Dirichlet == 2
		S.BC = 4; % Wire, Periodic in 1D, Dirichlet in 2D
	elseif n_Dirichlet == 3
		S.BC = 1; % Dirichlet in all 3D
	end
else
	% if user does not provide any BC, set default to periodic in 3D
	S.BC = 2;
	S.BCx = 0; S.BCy = 0; S.BCz = 0;
end 

L1 = S.L1; L2 = S.L2; L3 = S.L3;

% S.Nx, S.Ny, S.Nz is number of intervals now
if S.Nx > 0 && S.Ny > 0 && S.Nz > 0
	S.dx = S.L1 / S.Nx;
	S.dy = S.L2 / S.Ny;
	S.dz = S.L3 / S.Nz;
elseif S.ecut > 0
	S.mesh_spacing = Ecut2h(S.ecut, S.FDn);
	S.Nx = max(round(S.L1/S.mesh_spacing),1);
	S.Ny = max(round(S.L2/S.mesh_spacing),1);
	S.Nz = max(round(S.L3/S.mesh_spacing),1);
	S.dx = S.L1 / S.Nx;
	S.dy = S.L2 / S.Ny;
	S.dz = S.L3 / S.Nz;
elseif S.mesh_spacing > 0
	S.Nx = max(round(S.L1/S.mesh_spacing),1);
	S.Ny = max(round(S.L2/S.mesh_spacing),1);
	S.Nz = max(round(S.L3/S.mesh_spacing),1);
	S.dx = S.L1 / S.Nx;
	S.dy = S.L2 / S.Ny;
	S.dz = S.L3 / S.Nz;
end

dx = S.dx; dy = S.dy; dz = S.dz;
S.dV = S.dx * S.dy * S.dz * S.Jacb;

% Finite-difference discretization
S.Nx = S.Nx + S.BCx;
S.Ny = S.Ny + S.BCy;
S.Nz = S.Nz + S.BCz;
Nx = S.Nx; Ny = S.Ny; Nz = S.Nz;
S.N = S.Nx * S.Ny * S.Nz;

% map atom positions back to domain for periodic domain	
% in x direction	
isAtomOutx = sum(S.Atoms(:,1) < 0 | S.Atoms(:,1) >= S.L1) > 0;	
if (isAtomOutx)	
	if S.BCx == 0	
		S.Atoms(:,1) = mod(S.Atoms(:,1), S.L1);	
		for ityp = 1:S.n_typ	
			S.Atm(ityp).coords(:,1) = mod(S.Atm(ityp).coords(:,1),S.L1);	
		end	
		fprintf(' Warning: mapped atom position back to domain in x dir!\n');	
	else	
		error('Atom out of domain in x dir!');	
	end	
end	
% in y direction	
isAtomOuty = sum(S.Atoms(:,2) < 0 | S.Atoms(:,2) >= S.L2) > 0;	
if (isAtomOuty)	
	if S.BCy == 0	
		S.Atoms(:,2) = mod(S.Atoms(:,2), S.L2);	
		for ityp = 1:S.n_typ	
			S.Atm(ityp).coords(:,2) = mod(S.Atm(ityp).coords(:,2),S.L2);	
		end	
		fprintf(' Warning: mapped atom position back to domain in y dir!\n');	
	else	
		error('Atom out of domain in y dir!');	
	end	
end	
% in z direction	
isAtomOutz = sum(S.Atoms(:,3) < 0 | S.Atoms(:,3) >= S.L3) > 0;	
if (isAtomOutz)	
	if S.BCz == 0	
		S.Atoms(:,3) = mod(S.Atoms(:,3), S.L3);	
		for ityp = 1:S.n_typ	
			S.Atm(ityp).coords(:,3) = mod(S.Atm(ityp).coords(:,3),S.L3);	
		end	
		fprintf(' Warning: mapped atom position back to domain in z dir!\n');	
	else	
		error('Atom out of domain in z dir!');	
	end	
end	

if (isAtomOutx || isAtomOuty || isAtomOutz)	
	fprintf(' After mapping\n COORD:\n');	
	disp(S.Atoms);	
end

% Finite difference weights of the second derivative
FDn = S.FDn;
w2 = zeros(1,FDn+1) ;
for k=1:FDn
	w2(k+1) = (2*(-1)^(k+1))*(factorial(FDn)^2)/...
		(k*k*factorial(FDn-k)*factorial(FDn+k));
	w2(1) = w2(1)-2*(1/(k*k));
end
S.w2 = w2;

% Finite difference weights of the first derivative
w1 = zeros(1,FDn) ;
for k=1:FDn
	w1(k+1) = ((-1)^(k+1))*(factorial(FDn)^2)/...
		(k*factorial(FDn-k)*factorial(FDn+k));
end
S.w1 = w1;

% Weights for spatial integration over domain
% S.W = IntgWts(S.Nx,S.Ny,S.Nz,S.BCx,S.BCy,S.BCz,S.xin,S);
if S.cell_typ == 1 || S.cell_typ == 2
	S.W = ones(S.N,1) * (S.dx*S.dy*S.dz*S.Jacb);
else
	S.W = IntgWts(S.Nx,S.Ny,S.Nz,S.BCx,S.BCy,S.BCz,S.xin,S);
end

% assert(abs(sum(S.W)-pi*(S.xout*S.xout-S.xin*S.xin)*S.L3/25)<1e-6, ...
%      'Incorrect weights for spatial integration!');

% Brillouin-Zone Sampling
S = Generate_kpts(S);

% Create spherical harmonics for poisson solve for isolated clusters
if (S.BCx == 1 && S.BCy == 1 && S.BCz == 1)
	% Calculate Spherical Harmonics with origin shifted to the center of the domain
	xx_aug = (0-FDn:Nx+FDn-1)*dx;% - L1/2;
	yy_aug = (0-FDn:Ny+FDn-1)*dy;% - L2/2;
	zz_aug = (0-FDn:Nz+FDn-1)*dz;% - L3/2;
	[XX_AUG_3D,YY_AUG_3D,ZZ_AUG_3D] = ndgrid(xx_aug,yy_aug,zz_aug);
	% Find distances
	RR_AUG_3D = calculateDistance(XX_AUG_3D,YY_AUG_3D,ZZ_AUG_3D,L1/2,L2/2,L3/2,S);

	XX_AUG = reshape(XX_AUG_3D,[],1);
	YY_AUG = reshape(YY_AUG_3D,[],1);
	ZZ_AUG = reshape(ZZ_AUG_3D,[],1);
	RR_AUG = reshape(RR_AUG_3D,[],1);

	S.RR_AUG = RR_AUG;
	S.RR_AUG_3D = RR_AUG_3D;

	in_flag = ones(Nx+2*FDn,Ny+2*FDn,Nz+2*FDn);
	in_flag(1:FDn,:,:) = 0;
	in_flag(:,1:FDn,:) = 0;
	in_flag(:,:,1:FDn) = 0;
	in_flag(FDn+Nx+1:end,:,:) = 0;
	in_flag(:,FDn+Ny+1:end,:) = 0;
	in_flag(:,:,FDn+Nz+1:end) = 0;
	isIn = (in_flag~=0);

	S.isIn = isIn;
	S.RR = RR_AUG(isIn);

	pos_node_cart = coordinateTransformation(S,[XX_AUG,YY_AUG,ZZ_AUG],'noncart2cart_dis');
	pos_atm_cart = coordinateTransformation(S,[L1/2,L2/2,L3/2],'noncart2cart_dis');
	XX_AUG = bsxfun(@minus,pos_node_cart(:,1),pos_atm_cart(:,1));
	YY_AUG = bsxfun(@minus,pos_node_cart(:,2),pos_atm_cart(:,2));
	ZZ_AUG = bsxfun(@minus,pos_node_cart(:,3),pos_atm_cart(:,3));

	l_cut = 6;
	SH = repmat(struct([]),l_cut+1,1);
	for l = 0:l_cut
		for m = -l:l
			SH(l+1).Ylm_AUG(:,m+l+1) = sphericalHarmonics(XX_AUG,YY_AUG,ZZ_AUG,l,m,'real');
			Ylm_AUG_TEMP = SH(l+1).Ylm_AUG(:,m+l+1);
			SH(l+1).Ylm(:,m+l+1) = Ylm_AUG_TEMP(isIn);
			% SH(l+1).Ylm(:,m+l+1) = sphericalHarmonics(XX,YY,ZZ,l,m,'real');
		end
	end

	S.l_cut = l_cut;
	S.SH = SH;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set default values to the same default as SPARC %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first find effective mesh size
if S.cell_typ < 3
	dx2_inv = 1/(S.dx * S.dx);
	dy2_inv = 1/(S.dy * S.dy);
	dz2_inv = 1/(S.dz * S.dz);
	h_eff = sqrt(3.0 / (dx2_inv + dy2_inv + dz2_inv));
elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
	dx2_inv = 1/(S.dx * S.dx);
	dy2_inv = 1/(((S.xin+S.xout)/2)*S.dy)^2;
	dz2_inv = 1/(S.dz * S.dz);
	h_eff = sqrt(3.0 / (dx2_inv + dy2_inv + dz2_inv));
end
% find npl
if S.npl < 0
	fprintf('## Chebyshev polynomial degree not provided, finding npl ...\n');
	S.npl = Mesh2ChebDegree(h_eff);
	fprintf('## Based on the mesh size, npl is set to: %d\n',S.npl);
end

% Nev
if S.Nev < 0
	fprintf('## Number of states not provided, finding Nev ...\n');
	S.Nev = floor(S.Nelectron / 2) * 1.2 + 5; 
	S.Nev = round(S.Nev);
	fprintf('## Based on the number of electrons, Nev is set to: %d\n',S.Nev);
end

% SCF_tol
if S.SCF_tol < 0
	log10_neatom = log10(S.Nelectron / S.n_atm);
	target_force_accuracy = -1.0;
	target_energy_accuracy = -1.0;
	if S.accuracy_level >= 0
		target_force_accuracy = 10^(S.accuracy_level + 1);
	elseif S.target_force_accuracy > 0
		target_force_accuracy = S.target_force_accuracy;
	elseif S.target_energy_accuracy > 0
		target_energy_accuracy = S.target_energy_accuracy;
	end
	
	% if none of the accuracy levels are specified, set force_accuracy to
	% 1e-3
	if target_force_accuracy < 0  && target_energy_accuracy < 0 
		target_force_accuracy = 1e-3;
	end
	
	% choose SCF tol based on the desired accuracy
	if target_force_accuracy > 0
		a = 1.03;
		b = -0.33145 + 1.1761 * log10_neatom;
		S.SCF_tol = 10^((log10(target_force_accuracy) - b)/a);
	elseif target_energy_accuracy > 0
		a = 2.02;
		b = -0.25440 + 1.0000 * log10_neatom;
		S.SCF_tol = 10^((log10(target_energy_accuracy) - b)/a);
	end
	fprintf('## Based on the desired accuracy, SCF_tol is set to: %.3e\n',S.SCF_tol);
end

% poisson_tol
if S.poisson_tol < 0
	fprintf('## Poisson tolerance not provided, choosing poisson_tol ...\n')
	S.poisson_tol = S.SCF_tol * 0.01; 
	fprintf('## poisson_tol is set to: %.3e\n',S.poisson_tol);
end

% pseudocharge_tol
if S.pseudocharge_tol < 0
	fprintf('## Pseudocharge tolerance not provided, choosing pseudocharge_tol ...\n')
	S.pseudocharge_tol = S.SCF_tol * 0.01;
	fprintf('## pseudocharge_tol is set to: %.3e\n',S.pseudocharge_tol);
end

% % kerker_tol
% if S.kerker_tol < 0
%     S.kerker_tol = h_eff * h_eff * 0.001;
% end

% precond_tol
if S.precond_tol < 0
	S.precond_tol = h_eff * h_eff * 0.001;
end

% mixing parameter for simple mixing
if S.MixingParameterSimple < 0
	S.MixingParameterSimple = S.MixingParameter;
end

% Preconditioner for SCF convergence
if S.MixingVariable < 0
	% L_diag = sqrt(S.L1 * S.L1 + S.L2 * S.L2 + S.L3 * S.L3);
	% if (L_diag <= 20 || (S.BC ~= 2 && S.BC ~= 0))
	%     S.MixingVariable = 1;
	% else
	%     S.MixingVariable = 0;
	% end
	S.MixingVariable = 1; % set default mixing var to potential
end

if S.MixingPrecond < 0
	% if S.MixingVariable == 0
	%     S.MixingPrecond = 1;
	% else
	%     S.MixingPrecond = 0;
	% end
	S.MixingPrecond = 0; % set default precond to none
end

% set up coefficients
if S.MixingPrecond == 1 % kerker
	S.precondcoeff_a = 1.0;
	S.precondcoeff_lambda_TF = S.precond_kTF * S.precond_kTF;
	S.precondcoeff_k = 0;
elseif S.MixingPrecond == 2 % resta
	% put these in input options
	%max_q = sqrt(3) * 2*pi/0.08;  % Maximum q in the domain of fit
	max_q = 100;
	% mpower = 2; % The power used in the fit 
	% a0 = 0.25;   % parameters of the preconditioner
	% ktf = 1;     % parameters of the preconditioner
	% q0 = 1.36;  % parameters of the preconditioner
	% e0 = 5.7;   % parameters of the preconditioner
	% Rs = 2.76;  % parameters of the preconditioner
	mpower = S.precond_fitpow;
	ktf    = S.precond_kTF;
	a0     = S.precond_thresh;
	q0     = S.precond_resta_q0;
	Rs     = S.precond_resta_Rs;
	e0 = sinh(q0*Rs)/(q0*Rs);   % parameters of the preconditioner
	%[a_temp, lambda_TF_temp, const_temp] = fit_mixing_preconditioner(...
	%    100,1.36,5.7,2.76,0.25,1,2,2);
	[a_temp, lambda_TF_temp, const_temp] = fit_mixing_preconditioner(...
		max_q, q0, e0, Rs, a0, ktf, mpower, 2);
	% store coeffs in S
	S.precondcoeff_a         = a_temp;
	S.precondcoeff_lambda_TF = lambda_TF_temp;
	S.precondcoeff_k         = const_temp;
elseif S.MixingPrecond == 3 % truncated kerker
	% put these in input options
	%max_q = sqrt(3) * 2*pi/0.08;  % Maximum q in the domain of fit
	max_q = 100;
	%mpower = 2; % The power used in the fit 
	%a0 = 0.25;   % parameters of the preconditioner
	%ktf = 1;     % parameters of the preconditioner
	%q0 = 1.36;   % parameters of the preconditioner
	%Rs = 2.76;   % parameters of the preconditioner
	mpower = S.precond_fitpow;
	ktf    = S.precond_kTF;
	a0     = S.precond_thresh;
	q0     = S.precond_resta_q0;
	Rs     = S.precond_resta_Rs;
	e0 = sinh(q0*Rs)/(q0*Rs);   % parameters of the preconditioner
	[a_temp, lambda_TF_temp, const_temp] = fit_mixing_preconditioner(...
		max_q, q0, e0, Rs, a0, ktf, mpower, 1);
	S.precondcoeff_a         = a_temp;
	S.precondcoeff_lambda_TF = lambda_TF_temp;
	S.precondcoeff_k         = const_temp;
end

if S.nspin == 1
	S.occfac = 2;
else
	S.occfac = 1;
end

% Name of the relax file
% if S.RelaxFlag == 1
%     S.relaxfname = strcat(filename,'.geopt'); 
%     i = 1;
%     while exist(S.relaxfname,'file')
%         S.relaxfname = sprintf('%s.geopt_%d',filename,i);
%         i = i + 1;
%     end
% 
%     % if there are already 100 files, then start using .geopt only
%     OUT_MAX = 100;
%     if i > OUT_MAX
%         S.relaxfname = strcat(filename,'.geopt'); 
%     end
% end

% Name of the MD file
% if S.MDFlag == 1
%     S.mdfname = strcat(filename,'.aimd'); 
%     i = 1;
%     while exist(S.mdfname,'file')
%         S.mdfname = sprintf('%s.aimd_%d',filename,i);
%         i = i + 1;
%     end
% 
%     % if there are already 100 files, then start using .geopt only
%     OUT_MAX = 100;
%     if i > OUT_MAX
%         S.mdfname = strcat(filename,'.aimd'); 
%     end
% end


if (S.RelaxFlag || S.MDFlag)
	% Name of the restart file
	S.restartfname = strcat(filename,'.restart');
	% charge difference histories for charge extrapolation
	S.delta_rho_tm2 = zeros(S.N,1);
	S.delta_rho_tm1 = zeros(S.N,1);
	S.delta_rho_t   = zeros(S.N,1);
	S.delta_rho_in_tp1 = zeros(S.N,1);
	% S.dV_tm2 = S.dV;
	% S.dV_tm1 = S.dV;
	% S.dV_t   = S.dV;
	% S.dV_tp1 = S.dV;
	S.atom_pos_tm2  = S.Atoms_init;
	S.atom_pos_tm1  = S.Atoms_init;
	S.atom_pos_t    = S.Atoms_init;
	S.atom_pos_tp1  = S.Atoms_init;
	S.Atoms_old     = S.Atoms_init;
end

S.Relax_iter = 1;
S.ForceCount = 1;
S.amu2au = 1822.888485332371; % 1 au = 9.10938356e-31 Kg; 1 amu =  1.660539040e-27 Kg;
S.fs2atu = 41.34137333649300; %1atu = 2.418884326509e-17 s;

fprintf(' Creating differentiation matrices ...\n');
t1 = tic;

% Calculate discrete laplacian (1D) and discrete gradient indices' values
S = lapIndicesValues_1d(S);
S = gradIndicesValues(S);

% Calculate discrete laplacian
[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,[0 0 0]);
if S.cell_typ < 3
	S.Lap_std = S.lapc_T(1,1) * kron(speye(S.Nz),kron(speye(S.Ny),DL11))  +  S.lapc_T(2,2) * kron(speye(S.Nz),kron(DL22,speye(S.Nx))) + ...
				S.lapc_T(3,3) * kron(DL33,kron(speye(S.Ny),speye(S.Nx))) ;
	if (S.cell_typ == 2)
		MDL = S.lapc_T(1,2) * kron(speye(S.Nz),kron(DG2,DG1))  +  S.lapc_T(2,3) * kron(DG3,kron(DG2,speye(S.Nx))) + ...
			  S.lapc_T(1,3) * kron(DG3,kron(speye(S.Ny),DG1)) ;
		S.Lap_std = S.Lap_std + MDL;
	end
elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
	S.Lap_std = kron(speye(S.Nz),kron(speye(S.Ny),(DL11+DG1))) + kron(DL33,kron(speye(S.Ny),speye(S.Nx))) + ...
				kron(speye(S.Nz),kron(DL22,S.R2inv));
	if (S.cell_typ == 4 || S.cell_typ == 5)
		MDL = kron(speye(S.Nz),kron(DL22,speye(S.Nx))) +   kron(DG3,kron(DG2,speye(S.Nx)));
		S.Lap_std = S.Lap_std + MDL;
	end
end

% Calculate discrete gradient
S.grad_1 = blochGradient(S,[0 0 0],1);
S.grad_2 = blochGradient(S,[0 0 0],2);	
S.grad_3 = blochGradient(S,[0 0 0],3);

% Calculate preconditioners for negative discrete laplacian
[S.LapPreconL, S.LapPreconU] = ilu(S.Lap_std,struct('droptol',1e-5));

fprintf(' Done. (%.3f sec)\n', toc(t1));

% Estimate memory usage
S.memory_usage = estimate_memory(S);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function memory_usage = estimate_memory(S)
% estimate the memory required for the simulation
size_double = 8;

fd_order = S.FDn * 2;

N_ex = (S.Nx+fd_order) * (S.Ny+fd_order) * (S.Nz+fd_order);

% orbitals (dominant)
ncpy_orbitals = 4; % 4 copies required during chebyshev filtering
if S.nspin ~= 1, ncpy_orbitals = ncpy_orbitals * 2; end
% for kpoints, the factor 2 is for complex entries
if S.tnkpt ~= 1, ncpy_orbitals = ncpy_orbitals * 2 * S.tnkpt; end
memory_orbitals = S.N * S.Nev * size_double * ncpy_orbitals;

% sparse matrices
% Lap_std, grad_i, {GIi, GJi, GVi}
% (nnz * 2 + ncol + 1) * size_double
nnz_sparse = (2 * (3*fd_order+1) + 3 * fd_order) * S.N;
memory_sparse = (2*nnz_sparse + S.N + 1) * size_double ...
	+ 9*fd_order*S.N * size_double + (S.BC ~= 1) * 6*fd_order*S.N;

% vectors, rho, phi, Veff, ...; isIn, RR_AUG, RR_AUG_3D
memory_vectors = 14 * S.N * size_double + ...
	+ (S.BC == 1) * 3 * N_ex * size_double ...
	+ (S.RelaxFlag || S.MDFlag) * 4 * S.N * size_double;
    
% spherical harmonics
if S.BC == 1
	memory_SH = (6+1)^2 * (S.N + N_ex) * size_double;
else 
	memory_SH = 0;
end

% history matrix
memory_hist = 2 * S.MixingHistory * S.N * size_double;

% total
memory_usage = memory_orbitals + memory_sparse + memory_vectors ...
	+ memory_SH + memory_hist;

fprintf('\n');
fprintf(' Estimated memory usage:\n');
fprintf(' Total: %s\n', print_mem(memory_usage));
fprintf(' orbitals            : %s\n', print_mem(memory_orbitals));
fprintf(' sparse matrices     : %s\n', print_mem(memory_sparse));
fprintf(' global-size vectors : %s\n', print_mem(memory_vectors));
if S.BC == 1
	fprintf(' spherical harmonics : %s\n', print_mem(memory_SH));
end
fprintf(' mixing histories    : %s\n', print_mem(memory_hist));
fprintf('\n');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mem_str, nXB, XB] = print_mem(nBytes)
% given memory size in Bytes, print in appropriate units

if nBytes < 0
	error('Memory size must be non-negative!');
end
scale = floor(log(nBytes)/log(1024));
if nBytes == 0, scale = 0; end
switch scale
    case 0
		nXB = nBytes; XB = 'B';
    case 1
        nXB = nBytes/(1024); XB = 'kB';
    case 2
		nXB = nBytes/(1024^2); XB = 'MB';
    case 3
		nXB = nBytes/(1024^3); XB = 'GB';
    case 4
		nXB = nBytes/(1024^4); XB = 'TB';
    case 5
		nXB = nBytes/(1024^5); XB = 'PB';
	otherwise
		% use PB for all larger mem size
		nXB = nBytes/(1024^5); XB = 'PB';
end

if scale == 0
	mem_str = sprintf('%7.0f %s',nXB,XB);
else
	mem_str = sprintf('%7.2f %s',nXB,XB);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function npl = Mesh2ChebDegree(h) 
	% the relation between h and npl is fit with a cubic polynomial
	% p(x) = p3 * x^3 + p2 * x^2 + p1 * x + p0.
	p3 = -700. / 3.;
	p2 = 1240. / 3.;
	p1 = -773. / 3.;
	p0 = 1078. / 15.;
	if (h > 0.7) 
		npl = 14;
	else 
		npl = ((p3 * h + p2) * h + p1) * h + p0;
	end
	npl = round(npl);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = inpt_defaults()
% SET_DEFAULTS sets up the default parameters

% cell_typ: 1 - orthogonal unit cell or 2 - non-orthogonal unit cell
cell_typ = 1; % default is orthogonal

% lat_vec: Lattice unit vectors
lat_vec = eye(3);
% corresponding metric_T, grad_T, lapc_T and Jacb
metric_T = eye(3);
grad_T   = eye(3);
lapc_T   = eye(3);
Jacb     = 1;

% L1, L2, L3: domain size in each direction
% no defaults
L1 = 0.0;
L2 = 0.0;
L3 = 0.0;

% Nx, Ny, Nz: number of finite-difference grids in each direction
% no defaults
Nx = 0;
Ny = 0;
Nz = 0;
% N: N = Nx * Ny * Nz
N  = 0;

% Ecut
ecut = -1.0;

% mesh spacing
mesh_spacing = -1.0;

% kptgrid: k-point grid
kptgrid = [0.0 0.0 0.0]; % default is Gamma-point

% shift in k-point grid	
kptshift = [0.0 0.0 0.0];

% nkpt: number of k-points
nkpt = [1 1 1];

% tnkpt: number of Time-Reversal Symmetry reduced k-points
tnkpt = 1;

% wkpt: weights for k-points
wkpt = 1;

% BC: boundary conditions. 
% 1 - Zero Dirichlet, 2 - Periodic, 
% 3 - surface in x-y, 4 - wire along z-direction
BC = -1; % will be set up later after reading user inputs
BCx = -1; BCy = -1; BCz = -1;

% isBS: flag for band structure calculation
isBS = 0; % Default is off

% lattice: lattice type, for band structure calculation
lattice = 'undefined';

% SCF_tol: SCF tolerance
SCF_tol = -1.0; % default will be set after inpt is read

% Nev: Number of states/bands
Nev = -1; % default will be set after inpt is read

% poisson_tol: Poisson tolerance
poisson_tol = -1; % default will be set after inpt is read

% pseudocharge_tol: Pseudocharge (rb) tolerance
pseudocharge_tol = -1;

% Cst: Factor for conversion from Ha to eV
Cst = 27.21138602;

% Temp: Electronic temperature, beta := 1/ (kB * Temp)
%Temp = 315.7751307269723;
kB = (8.6173303e-5)/Cst;
%bet = 1.0 / (kB * Temp); % beta = 1 / smearing
elec_T_type = 1; % gaussian smearing
bet = -1;
Temp = -1;

% npl: Degree of Chebyshev polynomial
npl = -1; % default will be set after inpt is read

% conditioner: mixing type and preconditioner
% 1 - Potential mixing , 2 - Density + Kerker, 3 - Density mixing
% conditioner = 1;

% FDn: half of finite difference order
FDn = 6;

% rc_ref: rc reference
rc_ref = 0.5;

% max_relax_it: Maximun number of relaxations
max_relax_it = 0;

% dbg_switch: debug switch
dbg_switch = 0;

% xc: Exchange-correlation functional 
% 0 - LDA_PW, 1 - LDA_PZ, 2 - PBE(GGA)
xc = 0;

% Nelectron: number of electrons
% no default
Nelectron = 0;

% n_typ: number of atom types
n_typ = 0;

S = struct(...
	'cell_typ',cell_typ,'lat_vec',lat_vec,'metric_T',metric_T,...
	'grad_T',grad_T, 'lapc_T',lapc_T,'Jacb',Jacb,'L1',L1,'L2',L2,'L3',L3,...
	'Nx',Nx,'Ny',Ny,'Nz',Nz,'N',N,'ecut',ecut,'mesh_spacing',mesh_spacing,'kptgrid',kptgrid,...
	'kptshift',kptshift,'nkpt',nkpt,'tnkpt',tnkpt,'wkpt',wkpt,'BC',BC,'BCx',BCx,'BCy',BCy,'BCz',BCz,...
	'isBS',isBS,'lattice',lattice,'SCF_tol',SCF_tol,'Nev',Nev,'poisson_tol',poisson_tol,...
	'pseudocharge_tol',pseudocharge_tol, 'Cst',Cst,'kB',kB,'elec_T_type',elec_T_type,...
	'Temp',Temp,'bet',bet,'npl',npl,'FDn',FDn,...
	'rc_ref',rc_ref,'max_relax_it',max_relax_it,...
	'dbg_switch',dbg_switch,'xc',xc,...
	'Nelectron',Nelectron,'n_typ',n_typ);

S.TimeRevSym = 1; 

S.spin_typ = 0;
S.nspin = 1; % spin 

% EXTRA variables from SPARC
S.CheFSI_Optmz = 0;
S.chefsibound_flag = 0;
S.FixRandSeed = 0;
S.rhoTrigger = 4;
S.NetCharge = 0;
S.MAXIT_SCF = 100;
S.MINIT_SCF = 3;
S.MAXIT_POISSON = 1000;
S.accuracy_level = -1;
S.target_force_accuracy = -1.0;
S.target_energy_accuracy = -1.0;
S.TOL_RELAX = 5e-4;
S.TOL_LANCZOS = 1e-2;
% S.kerker_tol = -1;
S.precond_tol = -1;

S.precond_kTF = 1;
S.precond_thresh = 0.25;
S.precond_fitpow = 2;
S.precond_resta_q0 = 1.36;
S.precond_resta_Rs = 2.76;

% set back to -1
S.MixingVariable = 1;
S.MixingPrecond = 0;
S.Pf_guess = [];

S.MixingHistory = 7;
S.MixingParameter = 0.3;
S.MixingParameterSimple = -1.0; % for simple mixing, set up later
S.PulayFrequency = 1;
S.PulayRestartFreq = 0;
S.TWtime = 1000000000;
S.RelaxFlag = 0;
S.RelaxMeth = 'LBFGS';
S.max_relax_it = 100;
S.max_dilatation = 1.2;	
S.TOL_RELAX_CELL = 0.01; % in GPa (max pressure)
S.MDFlag = 0;
S.RestartFlag = 0;
S.MDMeth = 'NVE';
S.MD_dt = 1.0;
S.MD_Nstep = 0;
S.ion_T = -1.0;
S.thermos_TF = -1.0;
S.ion_elec_eqT = 1;
S.ion_vel_dstr = 2;
S.NLCG_sigma = 0.5;
S.qmass = 1;
S.L_history = 20;
S.L_finit_stp = 5e-3;
S.L_maxmov = 0.2;
S.L_autoscale = 1;
S.L_lineopt = 1;
S.L_icurv = 1.0;
S.FIRE_dt = 1.0;
S.FIRE_mass = 1.0;
S.FIRE_maxmov = 0.2;
S.Calc_stress = 0;
S.Calc_pres = 0;
S.PrintForceFlag = 1;         % flag for printing forces
S.PrintAtomPosFlag = 1;       % flag for printing atomic positions
S.PrintAtomVelFlag = 1;       % flag for printing atomic velocities
S.PrintElecDensFlag = 0;      % flag for printing final electron density
S.PrintEigenFlag = 0;         % Flag for printing final eigenvalues and occupations
S.PrintMDout = 1;             % Flag for printing MD output in a .aimd file
S.PrintRelaxout = 1;          % Flag for printing relax output in a .relax file
S.Printrestart = 1;           % Flag for printing output needed for restarting a simulation
S.Printrestart_fq = 1;        % Steps after which the output is written in the restart file

% Origin of the unit cell wrt some global origin
S.xin = 0;  
S.yin = 0;
S.zin = 0;

% Cychel
S.alph = 0.0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = Write_output_init(S, filename)

% open .out file
outfname = strcat(filename,'.out'); 
i = 1;
while exist(outfname,'file')
	outfname = sprintf('%s.out_%d',filename,i);
	i = i + 1;
end

suffixNum = i-1; % save suffix number, only used if suffixNum > 0

% if there are already 100 files, then start using .out only
OUT_MAX = 100;
if i > OUT_MAX
	outfname = strcat(filename,'.out'); 
	suffixNum = -1;
end

% create an output file and write initial variables
fileID = fopen(outfname,'w');
if (fileID == -1) 
	error('\n Cannot open file "%s"\n',outfname);
end 

start_time = fix(clock);
fprintf(fileID,'***************************************************************************\n');
fprintf(fileID,'*                      M-SPARC v1.0.0 (Oct 07, 2020)                      *\n');
fprintf(fileID,'*   Copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech   *\n');
fprintf(fileID,'*           Distributed under GNU General Public License 3 (GPL)          *\n');
fprintf(fileID,'*                Date: %s  Start time: %02d:%02d:%02d                  *\n',date,start_time(4),start_time(5),start_time(6));
fprintf(fileID,'***************************************************************************\n');
fprintf(fileID,'                           Input parameters                                \n');
fprintf(fileID,'***************************************************************************\n');

fprintf(fileID,'CELL: %f %f %f \n',S.L1,S.L2,S.L3);
if max(max(abs(S.lat_vec - eye(3)))) > 1e-14
	fprintf(fileID,'LATVEC:\n');
	fprintf(fileID,'%.15f %.15f %.15f \n',S.lat_vec(1,:));
	fprintf(fileID,'%.15f %.15f %.15f \n',S.lat_vec(2,:));
	fprintf(fileID,'%.15f %.15f %.15f \n',S.lat_vec(3,:));
end
fprintf(fileID,'FD_GRID: %d %d %d\n',S.Nx-S.BCx,S.Ny-S.BCy,S.Nz-S.BCz);
fprintf(fileID,'FD_ORDER: %d\n',S.FDn*2);
%fprintf(fileID,'BOUNDARY_CONDITION: %d\n',S.BC);
str_BC = ['P', 'D'];
fprintf(fileID,'BC:');
fprintf(fileID,' %s',str_BC(S.BCx+1));
fprintf(fileID,' %s',str_BC(S.BCy+1));
fprintf(fileID,' %s',str_BC(S.BCz+1));
fprintf(fileID,'\n');
if (S.BC==2 || S.BC==3 || S.BC==4)
	fprintf(fileID,'KPOINT_GRID: %d %d %d\n',S.nkpt);
	fprintf(fileID,'KPOINT_SHIFT: %d %d %d\n',S.kptshift);
end

if (S.spin_typ ~= 0) 
	fprintf(fileID,'SPIN_TYP: %d\n', S.spin_typ);  
end

if (S.elec_T_type == 0) 
	fprintf(fileID,'ELEC_TEMP_TYPE: fermi-dirac\n');  
elseif (S.elec_T_type == 1) 
	fprintf(fileID,'ELEC_TEMP_TYPE: gaussian\n');  
end
%fprintf(fileID,'ELEC_TEMP: %lf\n',S.elec_T);

fprintf(fileID,'SMEARING: %.9f\n',1/S.bet);
fprintf(fileID,'CHEB_DEGREE: %d\n',S.npl);
fprintf(fileID,'NSTATES: %d\n',S.Nev);
%fprintf(fileID,'NTYPES: %d\n',S.Ntypes);
fprintf(fileID,'EXCHANGE_CORRELATION: %s\n',S.XC);
fprintf(fileID,'CALC_STRESS: %d\n',S.Calc_stress);
if(S.Calc_stress == 0)
	fprintf(fileID,'CALC_PRES: %d\n',S.Calc_pres);
end
%if (S.MDFlag == 1 || S.RelaxFlag == 1)
%    fprintf(fileID,'TWTIME: %f\n',S.TWtime);
%end

if (S.CheFSI_Optmz == 1)
	fprintf(fileID,'CHEFSI_OPTMZ: %d\n',S.CheFSI_Optmz);
end
if (S.chefsibound_flag == 1)
	fprintf(fileID,'CHEFSI_BOUND_FLAG: %d\n',S.chefsibound_flag);
end
if (S.NetCharge ~= 0)
	fprintf(fileID,'NET_CHARGE: %d\n',S.NetCharge);
end
fprintf(fileID,'MAXIT_SCF: %d\n',S.MAXIT_SCF);
if (S.MDFlag == 1)
	fprintf(fileID,'MD_FLAG: %d\n',S.MDFlag);
	fprintf(fileID,'MD_METHOD: %s\n',S.MDMeth);
	fprintf(fileID,'MD_TIMESTEP: %.2f\n',S.MD_dt); 
	%fprintf(fileID,'ATOMIC_MASS:');
	%for (i = 0; i < S.Ntypes; i++)	 
	%     fprintf(fileID,' %.15f', S.Mass[i]);      
	%end
	fprintf(fileID,'MD_NSTEP: %d\n',S.MD_Nstep);
	fprintf(fileID,'ION_ELEC_EQT: %d\n',S.ion_elec_eqT);
	% fprintf(fileID,'ION_VEL_DSTR: %d\n',S.ion_vel_dstr);
	fprintf(fileID,'ION_TEMP: %f\n',S.ion_T);
	% if(strcmp(S.MDMeth,'NVT_NH'))
		% fprintf(fileID,'ION_TEMP_END: %lf\n',S.thermos_Tf);
		% fprintf(fileID,'QMASS: %lf\n',S.qmass);
	% end
end
if (S.RelaxFlag==1)
	fprintf(fileID,'RELAX_FLAG: %d\n',S.RelaxFlag);
	fprintf(fileID,'RELAX_METHOD: %s\n',S.RelaxMeth);
	fprintf(fileID,'RELAX_NITER: %d\n',S.max_relax_it);
	if(strcmp(S.RelaxMeth,'LBFGS'))
		fprintf(fileID,'L_HISTORY: %d\n',S.L_history);
		fprintf(fileID,'L_FINIT_STP: %f\n',S.L_finit_stp);
		fprintf(fileID,'L_MAXMOV: %f\n',S.L_maxmov);
		fprintf(fileID,'L_AUTOSCALE: %d\n',S.L_autoscale);
		fprintf(fileID,'L_LINEOPT: %d\n',S.L_lineopt);
		fprintf(fileID,'L_ICURV: %f\n',S.L_icurv);
	elseif (strcmp(S.RelaxMeth,'NLCG'))
		fprintf(fileID,'NLCG_SIGMA: %f\n',S.NLCG_sigma);
	elseif (strcmp(S.RelaxMeth,'FIRE'))
		fprintf(fileID,'FIRE_dt: %f\n',S.FIRE_dt);
		fprintf(fileID,'FIRE_mass: %f\n',S.FIRE_mass);
		fprintf(fileID,'FIRE_maxmov: %f\n',S.FIRE_maxmov);
	end
	fprintf(fileID,'TOL_RELAX: %.2E\n',S.TOL_RELAX);
elseif (S.RelaxFlag==2)	
	fprintf(fileID,'RELAX_FLAG: %d\n',S.RelaxFlag);	
	fprintf(fileID,'RELAX_NITER: %d\n',S.max_relax_it);	
	fprintf(fileID,'TOL_RELAX_CELL: %.2E\n',S.TOL_RELAX_CELL);	
	fprintf(fileID,'RELAX_MAXDILAT: %f\n',S.max_dilatation);    
end
fprintf(fileID,'TOL_SCF: %.2E\n',S.SCF_tol);
fprintf(fileID,'TOL_POISSON: %.2E\n',S.poisson_tol);
fprintf(fileID,'TOL_LANCZOS: %.2E\n',S.TOL_LANCZOS);
fprintf(fileID,'TOL_PSEUDOCHARGE: %.2E\n',S.pseudocharge_tol);
if (S.MixingVariable == 0)
	fprintf(fileID,'MIXING_VARIABLE: density\n');
elseif  (S.MixingVariable == 1)
	fprintf(fileID,'MIXING_VARIABLE: potential\n');
end
if (S.MixingPrecond == 0)
	fprintf(fileID,'MIXING_PRECOND: none\n');
elseif (S.MixingPrecond == 1)
	fprintf(fileID,'MIXING_PRECOND: kerker\n');
elseif (S.MixingPrecond == 2)
	fprintf(fileID,'MIXING_PRECOND: resta\n');
elseif (S.MixingPrecond == 3)
	fprintf(fileID,'MIXING_PRECOND: truncated_kerker\n');
end
if (S.MixingPrecond ~= 0)
	fprintf(fileID,'TOL_PRECOND: %.2E\n',S.precond_tol);
end
if (S.MixingPrecond == 1) % kerker
	%fprintf(fileID,'TOL_KERKER: %.2E\n',S.kerker_tol);
	fprintf(fileID,'PRECOND_KERKER_KTF: %.2f\n',S.precond_kTF);
elseif (S.MixingPrecond == 2) % resta
	%fprintf(fileID,'TOL_PRECOND: %.2E\n',S.precond_tol);
	fprintf(fileID,'PRECOND_RESTA_Q0: %.3f\n',S.precond_resta_q0);
	fprintf(fileID,'PRECOND_RESTA_RS: %.3f\n',S.precond_resta_Rs);
	fprintf(fileID,'PRECOND_FITPOW: %d\n',S.precond_fitpow);
elseif (S.MixingPrecond == 3) % truncated kerker
	fprintf(fileID,'PRECOND_KERKER_KTF: %.2f\n',S.precond_kTF);
	fprintf(fileID,'PRECOND_KERKER_THRESH: %.2f\n',S.precond_thresh);
	fprintf(fileID,'PRECOND_FITPOW: %d\n',S.precond_fitpow);
end

fprintf(fileID,'MIXING_PARAMETER: %.2f\n',S.MixingParameter);
if S.MixingParameter ~= S.MixingParameterSimple
	fprintf(fileID,'MIXING_PARAMETER_SIMPLE: %.2f\n',S.MixingParameterSimple);
end
fprintf(fileID,'MIXING_HISTORY: %d\n',S.MixingHistory);
fprintf(fileID,'PULAY_FREQUENCY: %d\n',S.PulayFrequency);
fprintf(fileID,'PULAY_RESTART_FREQ: %d\n',S.PulayRestartFreq);
fprintf(fileID,'REFERENCE_CUTOFF: %.2f\n',S.rc_ref);
fprintf(fileID,'RHO_TRIGGER: %d\n',S.rhoTrigger);
fprintf(fileID,'FIX_RAND: %d\n',S.FixRandSeed);
fprintf(fileID,'PRINT_FORCES: %d\n',S.PrintForceFlag);
fprintf(fileID,'PRINT_ATOMS: %d\n',S.PrintAtomPosFlag);
fprintf(fileID,'PRINT_EIGEN: %d\n',S.PrintEigenFlag);
fprintf(fileID,'PRINT_DENSITY: %d\n',S.PrintElecDensFlag);
if(S.MDFlag == 1)
	fprintf(fileID,'PRINT_MDOUT: %d\n',S.PrintMDout);
end
if(S.MDFlag == 1 || S.RelaxFlag == 1)  
	fprintf(fileID,'PRINT_VELS: %d\n',S.PrintAtomVelFlag);  
	fprintf(fileID,'PRINT_RESTART: %d\n',S.Printrestart);
	if(S.Printrestart == 1)
		fprintf(fileID,'PRINT_RESTART_FQ: %d\n',S.Printrestart_fq);
	end
end

if(S.RelaxFlag == 1)
	fprintf(fileID,'PRINT_RELAXOUT: %d\n',S.PrintRelaxout);
end

fprintf(fileID,'OUTPUT_FILE: %s\n',outfname);
if (S.RestartFlag == 1)
	fprintf(fileID,'RESTART_FLAG: %d\n',S.RestartFlag);
end

% fprintf(fileID,'***************************************************************************\n');
% fprintf(fileID,'                           Parallelization                                 \n');
% fprintf(fileID,'***************************************************************************\n');
% fprintf(fileID,'NP_KPOINT_PARAL: %d\n',S.npkpt);
% fprintf(fileID,'NP_BAND_PARAL: %d\n',S.npband);
% fprintf(fileID,'NP_DOMAIN_PARAL: %d %d %d\n',S.npNdx,S.npNdy,S.npNdz);
% fprintf(fileID,'NP_DOMAIN_PHI_PARAL: %d %d %d\n',S.npNdx_phi,S.npNdy_phi,S.npNdz_phi);

fprintf(fileID,'***************************************************************************\n');
fprintf(fileID,'                             Initialization                                \n');
fprintf(fileID,'***************************************************************************\n');
% fprintf(fileID,'Number of processors               :  %d\n',nproc);

if ( (abs(S.dx-S.dy) <=1e-12) && (abs(S.dx-S.dz) <=1e-12) ...
	&& (abs(S.dy-S.dz) <=1e-12) ) 
	fprintf(fileID,'Mesh spacing                       : % f (Bohr)\n',S.dx); 
else
	fprintf(fileID,'Mesh spacing in x-direction        : % f (Bohr)\n',S.dx); 
	fprintf(fileID,'Mesh spacing in y-direction        : % f (Bohr)\n',S.dy); 
	fprintf(fileID,'Mesh spacing in z direction        : % f (Bohr)\n',S.dz);     
end

if (S.BC==2 || S.BC==3 || S.BC==4) 
	fprintf(fileID,'Number of symmetry adapted k-points:  %d\n',S.tnkpt);       
end

fprintf(fileID,'Output printed to                  :  %s\n',outfname);

%if (S.PrintAtomPosFlag==1)
%    fprintf(fileID,'Atom positions printed to          :  %s\n',S.AtomFilename);      
%end

%if (S.PrintForceFlag==1)
%    fprintf(fileID,'Forces printed to                  :  %s\n',S.ForceFilename);
%end 

% if (S.PrintEigenFlag==1)
%     fprintf(fileID,'Final eigenvalues printed to       :  %s\n',S.EigenFilename);
% end

% if (S.MDFlag == 1 && S.PrintMDout == 1)
%     fprintf(fileID,'MD output printed to               :  %s\n',S.MDFilename);
% end

% if (S.RelaxFlag == 1 && S.PrintRelaxout == 1)
%     fprintf(fileID,'Relax output printed to            :  %s\n',S.RelaxFilename);
% end

fprintf(fileID,'Total number of atom types         :  %d\n',S.n_typ);
fprintf(fileID,'Total number of atoms              :  %d\n',S.n_atm);
fprintf(fileID,'Total number of electrons          :  %d\n',S.Nelectron);

for ityp = 1:S.n_typ
	fprintf(fileID,'Atom type %-2d (valence electrons)   :  %s %d\n',ityp,S.Atm(ityp).typ, S.Atm(ityp).Z);
	fprintf(fileID,'Pseudopotential                    :  %s\n',S.Atm(ityp).psdfname);     
	fprintf(fileID,'lloc                               :  %d\n',S.Atm(ityp).lloc);    
	if (S.MDFlag == 1)
		fprintf(fileID,'Atomic mass                        :  %.15f\n',S.Atm(ityp).Mass);
	end
	fprintf(fileID,'Pseudocharge radii of atom type %-2d :  %.2f %.2f %.2f\n',ityp,S.Atm(ityp).rb_x,S.Atm(ityp).rb_y,S.Atm(ityp).rb_z);
	fprintf(fileID,'Number of atoms of type %-2d         :  %d\n',ityp,S.Atm(ityp).n_atm_typ);
	% if (S.PrintAtomPosFlag == 1 && S.MDFlag == 0 && S.RelaxFlag == 0)
	%     fprintf(fileID,'Fractional coordinates of atoms of type %-2d    :\n',ityp);
	%     for j = 1:S.Atm(ityp).n_atm_typ
	%         fprintf(fileID,'%18.10f %18.10f %18.10f\n',S.Atm(ityp).coords(j,:)./[S.L1, S.L2, S.L3]);
	%     end
	% end
end

[~, mem_num, mem_unit] = print_mem(S.memory_usage);
fprintf(fileID, 'Estimated total memory usage       :  %-.2f %s\n', mem_num, mem_unit);

fclose(fileID);    



%-----------------------------------------
% Write atom positions to .static file 
%-----------------------------------------
if ((S.PrintAtomPosFlag == 1 || S.PrintForceFlag == 1) && S.MDFlag == 0 && S.RelaxFlag == 0)
	staticfname = strcat(filename,'.static'); 
	if suffixNum > 0
		staticfname = sprintf('%s.static_%d',filename,suffixNum);
	end
	% open file
	fid = fopen(staticfname,'w') ;
	assert(fid~=-1,'Error: Cannot open .static file %s',staticfname);

	if(S.PrintAtomPosFlag == 1)
		fprintf(fid,'***************************************************************************\n');
		fprintf(fid,'                            Atom positions                                 \n');
		fprintf(fid,'***************************************************************************\n');

		nFracCoord = sum(S.IsFrac);
		
		if nFracCoord == S.n_typ
			fprintf(fileID,'Fractional coordinates of atoms:\n');
			for j = 1:S.n_atm
				fprintf(fileID,'%18.10f %18.10f %18.10f\n',S.Atoms(j,:)./[S.L1, S.L2, S.L3]);
			end
		elseif nFracCoord == 0
			fprintf(fileID,'Cartesian coordinates of atoms (Bohr):\n');
			for j = 1:S.n_atm
				fprintf(fileID,'%18.10f %18.10f %18.10f\n',S.Atoms(j,:));
			end
		else
			for ityp = 1:S.n_typ
				if S.IsFrac(ityp)
					fprintf(fileID,'Fractional coordinates of %s:\n',S.Atm(ityp).typ); 
					for j = 1:S.Atm(ityp).n_atm_typ
						fprintf(fileID,'%18.10f %18.10f %18.10f\n',S.Atm(ityp).coords(j,:)./[S.L1, S.L2, S.L3]);
					end
				else
					fprintf(fileID,'Cartesian coordinates of %s (Bohr):\n',S.Atm(ityp).typ); 
					for j = 1:S.Atm(ityp).n_atm_typ
						fprintf(fileID,'%18.10f %18.10f %18.10f\n',S.Atm(ityp).coords(j,:));
					end
				end
			end
		end
	end

	% close file
	fclose(fid);
	
	S.staticfname = staticfname; % save to structure;
end


% save to structure
S.suffixNum   = suffixNum;
S.outfname    = outfname; 

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = Ecut2h(Ecut, FDn)
% Find mesh size h (Bohr) in finite-difference corresponding to a Ecut (Ha)
% in plane-wave codes.
% FDn: half of finite difference order

% this can be tuned to make the correspondence between h and Ecut more
% accurate
epsilon = 1e-1; % tolerance threshold for FD 2nd derivative approx.

% finite difference weights
w2 = zeros(1,FDn+1); 
for k=1:FDn
	w2(k+1) = (2*(-1)^(k+1))*(factorial(FDn)^2)/...
					(k*k*factorial(FDn-k)*factorial(FDn+k));
	w2(1) = w2(1)-2*(1/(k*k));
end

kk = linspace(0,pi,1001);
y_cos =  -w2(1) + (-2*w2(2:end)) * cos((1:FDn)' * kk);
freq_err = abs(y_cos - kk.^2);
kc = kk(find(freq_err < epsilon,1,'last'));

h = kc / sqrt(2*Ecut);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r4, p4, constant1] = fit_mixing_preconditioner(max_q, q0, e0, Rs, a0, ktf, m, choose_method)
% inputs (Defaults mentioned whenever available)
% max_q: max value of q on domain
% q0, e0, Rs: parmeters to Resta preconditioner
% m: power of polynomial used in the fit
% a0 = 0.25; % Threshold for Vasp preconditioner (AMIN in VASP code)
% ktf = 1;  % Thomas fermi length in Bohr
% choose_method = 1 for VASP-fit
%               = 2 for Resta-fit
% Use following values as default:
% q0 = 1.36, e0 = 5.7, Rs = 2.76 for diamond
% q0 = 1.10, e0 = 11.94, Rs = 4.28 for silicon
% q0 = 1.08, e0 = 16, Rs = 4.71 for germanium
%  values taken from https://journals.aps.org/prb/pdf/10.1103/PhysRevB.16.2717
Initial_guess_vasp{4} = [0.750005180027774  -0.612656341946256   0.176651619589811  -0.021091066193744   0.000902462058248  -0.000006309747662   0.520408099139247 -0.445424309904542   0.149745485915738  -0.061220832540094   0.013655050447764];
Initial_guess_vasp{3} = [ 0.749987658778237  -0.456416789370834   0.086223230682159  -0.005203450251595   0.000048194773012   0.716134813320535  -0.189056786409910 -0.071804972002278   0.039024958313563];
Initial_guess_vasp{2} = [ 0.750032410760875  -0.311559965943014   0.030528791849734  -0.000423605614679   0.938408853848147  -0.178908923603946   0.114299964666429];
Initial_guess_vasp{1} = [0.7499   -0.1849    0.0039    1.0273    0.3241];
Initial_guess_resta{4} = [0.824617865719867  -1.666939600095733   8.091267716160914  -7.570624248069551   1.103843271821802  -0.003945029200409   0.294295330483791 3.949649251298832  24.616714373390433 -30.280192648230070   4.491018917223143];
Initial_guess_resta{3} = [0.824427620366807  -1.323851625012163  -1.073719366292704   0.645483368849621   0.012403487988262   0.518438754779273  -2.408422023505795 -7.781456187967437   3.879228313224187];
Initial_guess_resta{2} = [0.824420567940696   0.390846073880153  -0.319378926237595  -0.006040724645454   2.592804423707606   2.985494672065439  -1.906590310047641];
Initial_guess_resta{1} = [0.8244    0.8065    0.0150    3.0937    4.5164];
N = 1000; % Number of Data points
q = (1e-10:max_q/N: max_q)';
if choose_method == 1
	precond_value = max(a0,(q.^2)./(q.^2 + ktf^2));   % expression for VASP preconditioner   
		f_rational =  fit(q.^2, (precond_value-a0), strcat('rat',num2str(m),num2str(m)),'StartPoint',Initial_guess_vasp{m-1}); % Fitting of rational polynomial with origin shift
elseif choose_method == 2
	precond_value = ((q0^2.*sin(q*Rs)./(e0*q*Rs) + q.^2))./(q0^2 + q.^2);   % expression for resta precnditioner
		f_rational = fit(q.^2, (precond_value-(1/e0)), strcat('rat',num2str(m),num2str(m)),'StartPoint',Initial_guess_resta{m-1}); % Fitting of rational polynomial with origin shift
end
coef_rational = coeffvalues(f_rational);  % coefficients of the fit

a = [1 coef_rational(1,m+2:end)];    % coefficents in denominator (increasing powers)
b = coef_rational(1,1:m) - coef_rational(m+1)*a(1,1:end-1)/a(1,end);    % coefficients in numerator (decreasing powers)
%[r,p,k1] = residue(b,a);     % partial fraction function
[r,p,~] = residue(b,a);     % partial fraction function

% Only slecting the real and one of the complex numbers in r and p 
r1 = r(imag(r)==0);
r2 = 2*r(imag(r)~=0);
r2 = sort(r2);
p1 = p(imag(p)==0);
p2 = p(imag(p)~=0);
p2 = sort(p2);

if isempty(r2) == 0
	p3 = zeros(floor(length(r2)/2),1);
	r3 = p3;
	for i = 1:length(r2)/2
		p3(i,1) = p2(2*i-1);
		r3(i,1) = r2(2*i-1);
	end
	p4 = [p1;p3];
	r4 = [r1;r3];
else
	p4 = p1;
	r4 = r1;
end
p4 = -p4; % here we make it consistent with kerker preconditioner

if choose_method ==1
	constant1 = a0+coef_rational(m+1)/a(1,end); % it is the constant which is added to the expansion
elseif choose_method == 2
	constant1 = 1/e0+coef_rational(m+1)/a(1,end); % it is the constant which is added to the expansion
end
%%% run this section only for printing the error
% if choose_method == 1
%     error_fit = max(abs(precond_value - f_rational(q.^2) - 1/e0));
% elseif choose_method == 2
%     error_fit = max(abs(precond_value - f_rational(q.^2) - a0));
% end
% display(error_fit)
%%% checking if the fit is correct
% part_f = zeros(length(q),1)+constant1;
% for i = 1:length(q)
%    for j = 1: length(r4)
%        if imag(r4(j)) ==0
%            part_f(i) = part_f(i)+r4(j)*((q(i))^2/((q(i))^2+p4(j)));
%        else
%             part_f(i) = part_f(i)+real(r4(j)*((q(i))^2/((q(i))^2+p4(j))));
%        end
%    end
% end
% plot(q, precond_value,'-*',q,part_f,'-*')
% xlim([0 10])
end

function [W] = IntgWts(Nx,Ny,Nz,BCx,BCy,BCz,xin,S)
	if S.cell_typ == 1 || S.cell_typ == 2
		W_x = ones(Nx,1)*S.dx;
		W_x(1) = W_x(1) * (1-BCx*0.5);
		W_x(Nx) = W_x(Nx) * (1-BCx*0.5);

		W_y = ones(Ny,1)*S.dy;
		W_y(1) = W_y(1) * (1-BCy*0.5);
		W_y(Ny) = W_y(Ny) * (1-BCy*0.5);

		W_z = ones(Nz,1)*S.dz;
		W_z(1) = W_z(1) * (1-BCz*0.5);
		W_z(Nz) = W_z(Nz) * (1-BCz*0.5);

		W = kron(W_z,kron(W_y,W_x)) * S.Jacb;
	elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
		W = IntgWts_cychel(Nx,Ny,Nz,BCx,BCy,BCz,xin,S);
	end
end


function [S] = Generate_kpts(S)
	nkpt = S.nkpt;
	if (S.BCx == 1 && nkpt(1) > 1)
		error(' nkpt cannot be greater than 1 in Dirichlet boundary direction (x)');
	end
	if (S.BCy == 1 && nkpt(2) > 1)
		error(' nkpt cannot be greater than 1 in Dirichlet boundary direction (y)');
	end
	if (S.BCz == 1 && nkpt(3) > 1)
		error(' nkpt cannot be greater than 1 in Dirichlet boundary direction (z)');
	end

	% Monkhorst-pack grid for Brillouin zone sampling
	MPG_typ1 = @(nkpt) (2*(1:nkpt) - nkpt - 1)/2; % MP grid points for infinite group order
	MPG_typ2 = @(nkpt) (0:nkpt-1); % MP grid points for finite group order

	if S.cell_typ < 3
		kptgrid_x = (1/nkpt(1)) * MPG_typ1(nkpt(1));
		kptgrid_y = (1/nkpt(2)) * MPG_typ1(nkpt(2));
		kptgrid_z = (1/nkpt(3)) * MPG_typ1(nkpt(3));
		sumx = 0;
		sumy = 0; 
		sumz = 0;
		% shift kpoint grid 
		kptgrid_x = kptgrid_x + S.kptshift(1) * (1/nkpt(1));
		kptgrid_y = kptgrid_y + S.kptshift(2) * (1/nkpt(2));
		kptgrid_z = kptgrid_z + S.kptshift(3) * (1/nkpt(3));
	
		% map k-points back to BZ
		temp_epsilon = eps; % include the right boundary k-points instead of left
		kptgrid_x = mod(kptgrid_x + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
		kptgrid_y = mod(kptgrid_y + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
		kptgrid_z = mod(kptgrid_z + 0.5 - temp_epsilon, 1) - 0.5 + temp_epsilon;
	elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
		kptgrid_x = (1/nkpt(1)) * MPG_typ1(nkpt(1));
		kptgrid_y = (1/nkpt(2)) * MPG_typ2(nkpt(2));
		kptgrid_z = (1/nkpt(3)) * MPG_typ1(nkpt(3));
		sumx = 0;
		sumy = nkpt(2); 
		sumz = 0;
	end    
	
	% Scale kpoints
	kptgrid_x = (2*pi/S.L1) * kptgrid_x;
	kptgrid_y = (2*pi/S.L2) * kptgrid_y;
	kptgrid_z = (2*pi/S.L3) * kptgrid_z;

	[kptgrid_X, kptgrid_Y, kptgrid_Z] = ndgrid(kptgrid_x,kptgrid_y,kptgrid_z);
	kptgrid = [reshape(kptgrid_X,[],1),reshape(kptgrid_Y,[],1),reshape(kptgrid_Z,[],1)];
	disp(' kpoint grid before symmetry:');
	disp(kptgrid);
	
	tnkpt = prod(nkpt);
	wkpt = ones(tnkpt,1)/tnkpt;% weights for k-points
	TOL = 1e-8;
	% Time-Reversal Symmetry to reduce k-points
	if S.TimeRevSym == 1
		Ikpt = zeros(tnkpt,1);
		Ikpt_rev = zeros(tnkpt,1);
		for ii = 1:tnkpt
			for jj = ii+1:tnkpt
				if (abs(kptgrid(ii,1) + kptgrid(jj,1) - sumx) < TOL) && (abs(kptgrid(ii,2) + kptgrid(jj,2) - sumy) < TOL) && (abs(kptgrid(ii,3) + kptgrid(jj,3) - sumz) < TOL)
					Ikpt(ii) = 1;
					Ikpt_rev(jj) = 1;
				end
			end
		end
		Ikpt = Ikpt>0.5;
		Ikpt_rev = Ikpt_rev>0.5;
		wkpt(Ikpt_rev) = 2*wkpt(Ikpt_rev);
		kptgrid = kptgrid(~Ikpt,:);
		wkpt = wkpt(~Ikpt);
		tnkpt = size(wkpt,1);
	end

	disp(' kpoint grid after symmetry:');	
	disp(kptgrid);
	% Store into the structure
	S.kptgrid = kptgrid;
	S.tnkpt   = tnkpt;
	S.wkpt    = wkpt;
end
