function out_struct = OneAtomSolver(input_struct, n_Uatm_typ, NLCC_flag, spin_Flag, XC, usefock, out_struct)
% @ brief     This code solves the 1D radial Kohn-Sham problem for a single atom.
%             1. Guess rho 
%             2. Solve Poisson Problem using rho.
%             3. Get XC from rho.
%             4. Solve Eigen Problem.
%             5. Update density (Periodic_Pulay Mixing)
%             6. GOTO 2 if density not converged.
%
% @ authors
%           Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%==========================================================================
format long 
t_tot = tic;

[filepath, ~, ~] = fileparts(which('msparc'));

fprintf("===========================================================\n")
fprintf("Launching SPARC-atom: One Atom Radial Kohn-Sham DFT solver\n")
fprintf("===========================================================\n")
%% Read Information from Pseudo-potential
% NOTE: Change the struct variable name 'S' to something else otherwise
% will interfere with the struct 'S' in the original M-SPARC

% NLCC
S.NLCC_flag = NLCC_flag;
S.tEigS = 0;
S.tSCF = 0;

%% Parameters
Nd = 400;                        % Number of discretization grid points
Nq = Nd;                         % Quadrature Grids
S.Nd = Nd; S.Nq = Nq;
S.alpha = 1;
S.beta = -0.45;
fprintf("Number of discretization grid points: \t%d\n", Nd)

Z = input_struct.Zatom;                 % Atomic Number
S.Zatom = Z;
Rmax = 20;                       % Expected half of domain (in Bohr), domain = [0,2Rmax]
S.MAXIT_SCF = 200;                 % The maximum number of scf cycles
S.SCF_tol = 1e-8;                  % SCF tolerance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choice of XC functionals available:

% 'LDA_PW'
% 'LDA_PZ' (only for spin - unpolarized)
% 'GGA_PBE'
% 'GGA_PBEsol'
% 'GGA_RPBE'
% 'SCAN'
% 'RSCAN'
% 'R2SCAN'
% 'PBE0' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(XC, 'R2SCAN')
    S.XC = 'RSCAN';
    fprintf("Changing XC from R2SCAN to RSCAN for convergence issues.\n");
elseif strcmp(XC, 'SCAN')
    S.XC = 'RSCAN';
    fprintf("Changing XC from SCAN to RSCAN for convergence issues.\n");
else
    S.XC = XC;
end

S.xc_rhotol = 1e-14;             % Density Tolerance for XC
S.xc_magtol = 1e-8;
S.xc_sigmatol = 1e-24;

% Spin-Polarized Flag
S.spinFlag = spin_Flag;

% SOC
S.nspinor = 1;
S.nspden = 1;

% hybrid functionals
S.usefock = usefock;
S.MAXIT_FOCK = -1;
S.MINIT_FOCK = -1;
S.FOCK_TOL = -1;
S.hyb_mixing = 0.25;
S.SCF_tol_init = -1;

% % check spin-orbit coupling
% if input_struct.pspsoc == 1
%     S.SOC_flag = 1;
% end

% spin un-polarized calculation
if ~S.spinFlag
    S.nspin = 1;
    S.nspden = 1;
    % if S.SOC_flag == 1
    %     S.nspinor = 2;
    % end
    fprintf("Spin un-polarized calculation.\n\n")
% collinear polarized calculation    
else
    S.nspin = 2;
    S.nspinor = 2;
    S.nspden = 2;
    fprintf("Spin polarized calculation.\n\n")
end

% Mixing Parameters
S.MixingHistory = 7;
S.MixingParameterSimple = 1;
S.MixingParameterSimpleMag = 1;
S.MixingParameter = 1;
S.MixingParameterMag = 1;
S.PulayRestartFlag = 0;

% Initialize neghrho flag
S.negrhoFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exchange correlation

%   Exchange:     "nox"    none                           iexch=0
%                 "slater" Slater (alpha=2/3)             iexch=1
%                 "pbex"   Perdew-Burke-Ernzenhof exch    iexch=2
%                       options: 1 -- PBE, 2 --PBEsol, 3 -- RPBE 4 --Zhang-Yang RPBE
%                 "rPW86x"  Refitted Perdew & Wang 86     iexch=3
%                 "scanx"  SCAN exchange                  iexch=4
%   
%   Correlation:  "noc"    none                           icorr=0
%                 "pz"     Perdew-Zunger                  icorr=1 
%                 "pw"     Perdew-Wang                    icorr=2
%                 "pbec"   Perdew-Burke-Ernzenhof corr    icorr=3
%                       options: 1 -- PBE, 2 --PBEsol, 3 -- RPBE
%                 "scanc"  SCAN correlation               icorr=4
%
%   Meta-GGA:     "nom"    none                           imeta=0
%                 "scan"   SCAN-Meta-GGA                  imeta=1
%                 "rscan"  rSCAN-Meta-GGA                 imeta=1
%                 "r2scan" r2SCAN-Meta-GGA                imeta=1
%
%   van der Waals "nov"    none                           ivdw=0
%                 "vdw1"   vdW-DF1                        ivdw=1
%                 "vdw2"   vdW-DF2                        ivdw=2

% decomposition of XC, ixc = [iexch,icorr imeta ivdw]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(S.XC, 'LDA_PW')
    S.xc = 0;
    S.isGradient = 0;
    S.ixc = [1 2 0 0];
elseif strcmp(S.XC, 'LDA_PZ')
    S.xc = 1;
    S.isGradient = 0;
    S.ixc = [1 1 0 0];
elseif strcmp(S.XC, 'GGA_PBE')
    S.xc = 2;
    S.isGradient = 1;
    S.ixc = [2 3 0 0];
    S.xc_option = [1 1];
elseif strcmp(S.XC, 'GGA_PBEsol')
    S.xc = 2;
    S.isGradient = 1;
    S.ixc = [2 3 0 0];
    S.xc_option = [2 2];
elseif strcmp(S.XC, 'GGA_RPBE')
    S.xc = 2;
    S.isGradient = 1;
    S.ixc = [2 3 0 0];
    S.xc_option = [3 3];
elseif strcmp(S.XC, 'SCAN') 
    S.xc = 4;
    S.SCF_tol = 1e-4;
    S.isGradient = 1;
    S.ixc = [4 4 1 0];
elseif strcmp(S.XC, 'RSCAN')
    S.xc = 4;
    S.SCF_tol = 1e-6;
	S.ixc = [5 5 1 0]; % 5: rscanx; 5: rscanc; 1: need kinetic energy density; 0: no vdWDF
    S.isGradient = 1; 
elseif strcmp(S.XC, 'R2SCAN') 
    S.xc = 4;
    S.SCF_tol = 1e-4;
    S.isGradient = 1;
    S.ixc = [6 6 1 0];
elseif strcmp(S.XC, 'PBE0') 
    S.xc = 41;
    S.usefock = 1;
    S.ixc = [2 3 0 0];
    S.xc_option = [1 1];
    S.isGradient = 1;
end

if (S.ixc(3) == 1 && S.NLCC_flag == 1)
		error('ERROR: Currently metaGGA functionals (SCAN, RSCAN, R2SCAN) do not support nonlinear core correction pseudopotential.');
elseif S.ixc(3) == 1
    if ispc % Windows system
        addpath(fullfile(filepath,'atom\mgga_atom'));
        % addpath('mgga_atom\');
    else % Mac/Linux
        addpath(fullfile(filepath,'atom/mgga_atom'));
        % addpath('mgga_atom/');
    end
elseif S.usefock == 1
    if ispc % Windows system
        addpath('exx_atom\');
    else % Mac/Linux
        addpath('exx_atom/');
    end
    S = exx_initialization_atom(S);
end

fprintf("Atomic Number of element: \t %d\n",Z)
fprintf("Exchange-Correlation used: \t %s\n",S.XC)
%% Basic Setup
% Exponential Grid
alpha = S.alpha;
beta = S.beta;

xmax = (alpha - alpha*exp(beta*2*Rmax))/2;
[D, x] = chebD(Nd,xmax);         % First derivative, grid r
w = xmax*clencurt(Nq);           % Clenshaw-Curtis weights of integration (in x)
r = (1/beta)*log(1 - x./alpha);  % r = r(x)
r(1) = 2*Rmax;
int_scale = -beta*(alpha - x);   % dr = int_scale*dx

% Laplacian matrix
L = (int_scale.^2).*D^2 + beta*int_scale.*D;     % d2/dr2  
S.Laplacian.matrix = L;          % On grid points in x

% Gradient matrix
grad_mat = int_scale.*D;            % d/dr
S.Gradient.matrix = grad_mat;

S.r = r;
S.x = x;
S.xmax = xmax;
S.w = w;
S.int_scale = int_scale;

% Guess rho from the pseudo-potential
r_grid_rho = input_struct.r_grid_rho;
rho_psp = input_struct.rho_isolated_guess;
r_rho_cut = max(r_grid_rho);

% Extrapolate guess rho using last two points
guessN = rho_psp(end);
guessN_1 = rho_psp(end-1);
rN = r_grid_rho(end);
rN_1 = r_grid_rho(end-1);
C1 = log(guessN/guessN_1)/(rN - rN_1);
C2 = guessN/exp(C1*rN);

index = r < r_rho_cut;
r_rho_in = r(index);
rho_guess = zeros(Nd+1,1);
rho_guess(index) = interp1(r_grid_rho, rho_psp, r_rho_in, 'spline');
rho_correct = zeros(Nd+1,1);
rho_correct(~index) = C2*exp(C1*r(~index));
rho_guess = rho_guess + rho_correct;
S.rho = zeros(Nd-1,3);
S.rho(:,1) = rho_guess(2:Nd);


% NLCC
rTilde = input_struct.r_grid_rho_Tilde;
rho_Tilde_psp = input_struct.rho_Tilde;
rTilde_cut = max(rTilde);
index = r < rTilde_cut;
rTilde_in = r(index);
rho_Tilde = zeros(Nd+1,1);
rho_Tilde(index) = interp1(rTilde, rho_Tilde_psp, rTilde_in, 'spline');
S.rho_Tilde = rho_Tilde(2:Nd);

% Account atomic states
Occ = getAtomicStates(Z);
n = Occ(1,:);
l = Occ(2,:);
f_tot = Occ(3,:)+Occ(4,:);
f_up = Occ(3,:); f_dw = Occ(4,:);
S.AtmStates.n = n;
S.AtmStates.l = l;
S.AtmStates.f_tot = f_tot;
S.AtmStates.f_up = f_up; S.AtmStates.f_dw = f_dw;

% Rho_up and Rho_dw
net_mag = abs(sum(f_up) - sum(f_dw));
integrand = [0;(S.r(2:S.Nd).^2).*S.rho(:,1);0];
scal = 4*pi*S.w*(integrand./S.int_scale);
if Z == 1
    S.rho(:,2) = S.rho(:,1);
    S.mag = S.rho(:,2) - S.rho(:,3);
else
    frac = (1+net_mag/scal)/2;
    S.rho(:,2) = frac*S.rho(:,1);
    S.rho(:,3) = (1-frac)*S.rho(:,1);
    S.mag = S.rho(:,2) - S.rho(:,3);
end
rho_guess_up = S.rho(:,2);
rho_guess_dw = S.rho(:,3);


%% Calculate Non-Local Projector Matrix (V_nl)
tnon_loc = tic;
lmax = input_struct.lmax;
nl_info = input_struct.Pot;
rc = input_struct.rc;
r_grid_vloc = input_struct.r_grid_vloc;

for i = 1:lmax+1
    gamma_Jl = nl_info(i).gamma_Jl;
    proj = nl_info(i).proj;       
    rc = input_struct.rc_max(i);
    index = r <= rc;
    r_proj_in = r(index);    
    V = zeros(Nd+1,Nd+1);
    for pp = 1 : input_struct.nproj(i)
        Chi_l = zeros(Nd+1,1);
        Chi_l(index) = interp1(r_grid_vloc,proj(:,pp),r_proj_in,'spline');
        rChi_l = (Chi_l.*r);
        wt_rChi_l = zeros(Nd+1,1);
        wt_rChi_l(index) = (w(index)'.*rChi_l(index))./int_scale(index);        
        V = V + gamma_Jl(pp)*(rChi_l*wt_rChi_l');        
    end    
    S.Vnl(i).matrix = V(2:Nd,2:Nd);    
end

tnon_loc = toc(tnon_loc);
fprintf("Time for storing V_nl:\t\t %0.10f s\n", tnon_loc)

%% Calculate the local part of the pseudo-potential
t_loc = tic;
% VJ from pseudo-potential
rcut = input_struct.r_grid_vloc(end);   % Cut-off radius for pseudo-potential
index = r < rcut;                % Stores logical values where r < rcut
r_in = r(index);                 % Stores the r<rcut values
r_grid_vloc = input_struct.r_grid_vloc;
Vloc = input_struct.Vloc;

% Spline interpolation of VJ to spectral grid
VJ = zeros(Nd+1,1);
VJ(~index) = -input_struct.Z;
VJ(index) = interp1(r_grid_vloc, r_grid_vloc.*Vloc, r_in, 'spline');
d2VJ = L*VJ;
VJ = VJ(2:Nd);
VJ = VJ./r(2:Nd);
VJ = diag(VJ);
S.d2VJ = d2VJ(2:Nd);
t_loc = toc(t_loc);
fprintf("Time for storing V_ext:\t\t %0.10f s\n", t_loc)

S.VJ.matrix = VJ;

%% Main Code
t_start = tic;
S = scf_loop(S, input_struct);
t_SCF = toc(t_start);S.tSCF = t_SCF;
fprintf("Time for SCF Cycle \t\t: %0.10f sec\n", t_SCF)

% Exact exchange potential 
if S.usefock == 1
    S.usefock = S.usefock+1;
end

% Exact exchange potential parameters
denMat_prev = densityMatrix(S);
count_Exx = 1;

% Exact Exchange Outer Loop
while(count_Exx <= S.MAXIT_FOCK)
    fprintf("*************************************************************\n")
    fprintf("<strong>No.%d Exx outer loop</strong>.\n",count_Exx)
    fprintf("*************************************************************\n")
    
    % Store orbitals and occupations for outer loop
    S.SCF_orbitals.matrix = S.orbitals.matrix;
    S = exx_operator(S);

    t_start = tic;
    S = scf_loop(S, input_struct);
    t_SCF = toc(t_start); 
    fprintf("Time for SCF Cycle \t\t: %0.10f sec\n", t_SCF)
    
    denMat = densityMatrix(S);
    err_Exx = abs(norm(denMat) - norm(denMat_prev))/norm(denMat_prev);
    fprintf('\nExx outer loop error\t\t: %.14e \n',err_Exx) ;
    if err_Exx < S.FOCK_TOL && count_Exx >= S.MINIT_FOCK
        break;
    else
        denMat_prev = denMat;
    end

    count_Exx = count_Exx + 1;
end

if S.usefock > 1
    if count_Exx > S.MAXIT_FOCK && err_Exx > S.FOCK_TOL
        disp('Warning: Exact Exchange outer loop did not converge. Maximum iterations reached!');
    else
        fprintf("Density Matrix has already converged to %0.4e.\n",S.FOCK_TOL)
        fprintf('<strong>Finished outer loop in %d steps!\n</strong>', count_Exx);
    end
end


% Display Energy
if S.usefock > 1
    S.SCF_orbitals.matrix = S.orbitals.matrix;
    S = evaluateExxEnergy(S);
end
S = EvaluateEnergyAtom(S);

n_final = [(S.n0)';(S.n1)';(S.n2)';(S.n3)'];
l_final = [zeros(size((S.n0)'));...
    ones(size((S.n1)'));...
    2*ones(size((S.n2)'));...
    3*ones(size((S.n3)'))];

% Sort
eigen_val_up = S.EigVal(:,1);
eigen_val_dw = S.EigVal(:,2);
[~,i] = sort(n_final);
n_final = n_final(i);
l_final = l_final(i);
eigen_val_up = eigen_val_up(i);
eigen_val_dw = eigen_val_dw(i);

fprintf("============================================================\n")
fprintf("\t\t <strong>Eigenvalues (Hartree)</strong>\n")
fprintf("=============================================================\n")
fprintf("<strong>n\t l\t   s\t\tOcc\t   Eigenvalue</strong>\n")
fprintf("=============================================================\n")

for i = 1:length(n_final)
    fprintf("%d\t %d\t +%0.1f\t\t %d\t %2.15f\n",n_final(i),l_final(i),0.5,S.f_up(i),eigen_val_up(i));
    fprintf("%d\t %d\t -%0.1f\t\t %d\t %2.15f\n",n_final(i),l_final(i),0.5,S.f_dw(i),eigen_val_dw(i));
end

t_tot = toc(t_tot);
fprintf("*************************************************************\n")
fprintf("\t\t\tTiming Info\n")
fprintf("*************************************************************\n")
fprintf('Total walltime\t\t\t: %0.10f sec\n',t_tot);
fprintf("_____________________________________________________________\n")

% Store the orbital-info to be used in main-code
out_struct(n_Uatm_typ).nmax_l = [length(S.n0) length(S.n1) length(S.n2) length(S.n3)];
out_struct(n_Uatm_typ).lmax = max(S.orbital_l.matrix(:));
out_struct(n_Uatm_typ).r = flip(S.r(2:S.Nd)); 
out_struct(n_Uatm_typ).rho = flip(S.rho);

if S.spinFlag
    store_Occ = [S.f0_up S.f1_up S.f2_up S.f3_up S.f0_dw S.f1_dw S.f2_dw S.f3_dw];
else
    store_Occ = [S.f0_tot S.f1_tot S.f2_tot S.f3_tot];
end

Nstates_half = size(S.occ.matrix,1);
for spinor = 1 : S.nspinor
    shift = (spinor - 1)*Nstates_half;
    states = 1;
    while(states <= Nstates_half)
        store_shift = 0;
        for l = 0 : out_struct(n_Uatm_typ).lmax
            nmax = out_struct(n_Uatm_typ).nmax_l(l+1);
            shift_proj = (spinor - 1)*nmax;
            % out_struct(n_Uatm_typ).orb(l+1).U = out_struct(n_Uatm_typ).Uval;
            for l_count = 1:nmax
                out_struct(n_Uatm_typ).orb(l+1).proj(:,l_count+shift_proj) = flip(S.orbitals.matrix(:,states+shift)./S.r(2:S.Nd));
                out_struct(n_Uatm_typ).orb(l+1).count(l_count+shift_proj) = store_Occ(store_shift+l_count+shift_proj);
                states = states + 1;
            end
            store_shift = store_shift + nmax;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Poisson Solve
function phi = poissonSolve(S,input_struct)
% Inputs:
% r = grid nodes
% D = Chebyshev differentiation Matrix
% rho = electron density(except at first point)

% Output:
% phi = electrostatic potential

Nd = S.Nd;
r = S.r;

% Left Hand Side
L = S.Laplacian.matrix(1:Nd,1:Nd);
L(1,:) = zeros(1,Nd);
L(1,1) = 1;

rho = S.rho(:,1);
rho = [0; rho];
r_rho = -4*pi*r(1:end-1).*rho;
RHS = [input_struct.Z; r_rho(2:end)];
rphi = L\RHS;
phi = rphi./r(1:Nd);
phi = phi(2:end);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Eigen Solve
function [S, eigvec, eigval] = eigenSolve(phi, Vxc, l, S)
% Inputs:
% r = grid nodes
% D = Chebyshev differentiation Matrix
% phi = electrostatic potential (except at first and last grid point)
% Vxc = Exchange Correlation Potential (except at first and last grid
% point)
% l = azimuthal quantum number

% Outputs:
% eigval = vector of eigen values
% eigvec = matrix containing all the eigen vectors

Nd = S.Nd;
r = S.r(2:Nd);

Vl = 0.5*l*(l+1)./(r.^2);
Vl = diag(Vl);
phi = diag(phi);

if ~S.spinFlag
    Vxc_up = diag(Vxc);
else
    Vxc_up = Vxc(:,1);
    Vxc_up = diag(Vxc_up);
    Vxc_dw = Vxc(:,2);
    Vxc_dw = diag(Vxc_dw);
end

% Hamiltonian
H_up = -0.5*S.Laplacian.matrix(2:Nd,2:Nd) + phi + Vxc_up + S.VJ.matrix ...
    + Vl + S.Vnl(l+1).matrix;

if S.spinFlag
    H_dw = -0.5*S.Laplacian.matrix(2:Nd,2:Nd) + phi + Vxc_dw + S.VJ.matrix ...
        + Vl + S.Vnl(l+1).matrix;
end

if S.xc == 4 && S.countPotential > 0 % mGGA    
    spin = 0.5;
    H_upx = @(x) H_up*x + evaluateMggaPotential_atom(S, l, x, spin);
elseif S.xc == 41 && S.usefock > 1 % Hybrid
    H_up = H_up + S.hyb_mixing*S.l_channel(l+1).Vexx_up;
end

% Up-spin
if S.xc == 4 && S.countPotential > 0 % mGGA
    [eigvec_up, eigval_up] = eigs(H_upx, Nd - 1, Nd - 1, 'sr',...
        'Tolerance',1e-14);
else
    [eigvec_up, eigval_up] = eig(H_up);
end

eigval_up = diag(eigval_up);

% Sorting the eigenvalues and eigenvectors
[~,i] = sort(eigval_up);
eigval_up = eigval_up(i);
eigvec_up = eigvec_up(:,i);

% Dw-spin
if ~S.spinFlag
    eigvec = [eigvec_up eigvec_up];
    eigval = [eigval_up eigval_up];
else
    
    if S.xc == 4 && S.countPotential > 0 % mGGA
        spin = -0.5;
        H_dwx = @(x) H_dw*x + evaluateMggaPotential_atom(S, l, x, spin);
    elseif S.xc == 41 && S.usefock > 1 % Hybrid
        H_dw = H_dw + S.hyb_mixing*S.l_channel(l+1).Vexx_dw;
    end
    
    if S.xc == 4 && S.countPotential > 0 % mGGA
        [eigvec_dw, eigval_dw] = eigs(H_dwx, Nd - 1, Nd - 1, 'sr',...
            'Tolerance',1e-14);
    else
        [eigvec_dw, eigval_dw] = eig(H_dw);
    end
    eigval_dw = diag(eigval_dw);
    
    % Sorting the eigenvalues and eigenvectors
    [~,i] = sort(eigval_dw);
    eigval_dw = eigval_dw(i);
    eigvec_dw = eigvec_dw(:,i);
    eigvec = [eigvec_up eigvec_dw];
    eigval = [eigval_up eigval_dw];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Anderson Weigted Average
function [x_wavg, f_wavg] = andersonWtdAvg(x_k, f_k, X, F)
% @brief    ANDERSONWTDAVG finds the weighted averages
%           x_wavg := x_k - X*Gamma
%           f_wavg := f_k - F*Gamma
%           where Gamma = inv(F' * F) * F' * f_k, F' denotes the transpose 
%           of F.
%
% @param x_k   Current guess/input vector.
% @param f_k   Current preconditioned residual.
% @param X     Iterate histories, X = [x_{k-m+1}-x_{k-m}, ..., x_k-x_{k-1}].
% @param F     Residual histories, F = [f_{k-m+1}-f_{k-m}, ..., f_k-f_{k-1}].
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech

Gamma = pinv(F'*F)*(F'*f_k);  
x_wavg = x_k - X * Gamma;
f_wavg = f_k - F * Gamma;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Periodic Pulay mixing
function [S,x_kp1] = Periodic_Pulay(S,g_k,x_k,iter,input_struct)
% For spin-unpolarized

% @brief    MIXING performs mixing of the previous SCF iterates to 
%           accelerate SCF convergence.
%
% @param g_k    Current output mixing iterate.
% @param x_k    Current input mixing iterate.
% @param iter   Iteration number.

m = S.MixingHistory;                        % Mixing History
p = 1;                                      % Pulay Frequency
beta = S.MixingParameter;                   % Mixing Parameter
omega = S.MixingParameterSimple;            % Simple mixing

beta_mag = S.MixingParameterMag;            % Magnetization Mixing Parameter
omega_mag = S.MixingParameterSimpleMag;     % Magnetization Simple Mixing

Pulay_mixing_flag = (rem(iter,p) == 0 && iter > 1);

if Pulay_mixing_flag
    amix = beta; amix_mag = beta_mag;
else
    amix = omega; amix_mag = omega_mag;
end

f_k = g_k - x_k;
if iter > 1
	f_km1 = S.mixing_hist_fkm1;
	x_km1 = S.mixing_hist_xkm1;
end

% store residual & iteration history
if iter > 1
	i_hist = mod(iter-2,m)+1;
	if (S.PulayRestartFlag ~= 0 && i_hist == 1)
		S.X = zeros(size(S.X)); S.F = zeros(size(S.F));
		S.X(:,1) = x_k - x_km1;
		S.F(:,1) = f_k - f_km1;
	else
		S.X(:,i_hist) = x_k - x_km1;
		S.F(:,i_hist) = f_k - f_km1;
	end
end

% apply Anderson extrapolation every p iters
if Pulay_mixing_flag
	% find weighted averages x_wavg, f_wavg
	[x_wavg, f_wavg] = andersonWtdAvg(x_k, f_k, S.X, S.F);
else 
	% simple mixing
	x_wavg = x_k; f_wavg = f_k;
end

if ~S.spinFlag
    f_tot = f_wavg; % for spin-unpolarized calculations, f_tot is just f_wavg
else
    f_tot = f_wavg(1:S.Nd - 1) + f_wavg(S.Nd : end);
    f_mag = f_wavg(1:S.Nd - 1) - f_wavg(S.Nd : end);
    sum_f_tot = sum(f_tot);
    sum_f_mag = sum(f_mag);
end

Pf = amix * f_tot; % No pre-conditioner

if S.spinFlag
    Pf(:,2) = amix_mag * f_mag; % No pre-conditioner
end

if S.spinFlag 
    sum_Pf_tot = sum(Pf(:,1));
    shift_Pf_tot = (sum_f_tot - sum_Pf_tot)/(S.Nd-1);
    Pf(:,1) = Pf(:,1) + shift_Pf_tot;
    
    sum_Pf_mag = sum(Pf(:,2));
    shift_Pf_mag = (sum_f_mag - sum_Pf_mag)/(S.Nd-1);
    Pf(:,2) = Pf(:,2) + shift_Pf_mag;
end

if ~S.spinFlag
    x_kp1 = x_wavg + Pf;
else
    x_kp1 = x_wavg + vertcat(Pf(:,1) + Pf(:,2), Pf(:,1) - Pf(:,2))/2;
end

negrho_count = sum(x_kp1 < 0);
if negrho_count > 0
    fprintf('Density got negative\n\n');
    S.negrhoFlag = 1;
end
x_kp1(x_kp1 < 0) = 0;

if negrho_count > 0
    func = [0;(S.r(2:S.Nd).^2).*x_kp1(1:S.Nd-1);0];
    if S.spinFlag
        func = [func./S.int_scale; func./S.int_scale];
    else
        func = func./S.int_scale;
    end
    scal = 4*pi*sum(S.w*reshape(func,S.Nd+1,S.nspin),2);
    x_kp1 = (input_struct.Z/scal)*x_kp1;
end

% update the history vectors
S.mixing_hist_fkm1 = f_k;
S.mixing_hist_xkm1 = x_k;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCF Loop
function S = scf_loop(S,input_struct)
%--------------------------------------------------------------------------
Nd = S.Nd;
Nq = S.Nq;
w = S.w;
r = S.r;
x = S.x;
xmax = S.xmax;
int_scale = S.int_scale;

n = S.AtmStates.n;
l = S.AtmStates.l;
f_tot = S.AtmStates.f_tot;
f_up = S.AtmStates.f_up; f_dw = S.AtmStates.f_dw;

if S.usefock == 1
    scf_tol = S.SCF_tol_init;
else
    scf_tol = S.SCF_tol;
end

%--------------------------------------------------------------------------
% Core occupations removed
Z_core = input_struct.Zatom - input_struct.Z;
Z_sum = 0;

for i = 1:length(f_tot)
    if Z_sum < Z_core
        Z_sum = Z_sum + f_tot(i);
        f_tot(i) = 100;
        f_up(i) = 100; f_dw(i) = 100;
        l(i) = 100;
        n(i) = 100;
    else
        break
    end
end
    
n = n(n~=100);
f_tot = f_tot(f_tot~=100);
f_up = f_up(f_up~=100); f_dw = f_dw(f_dw~=100);
l = l(l~=100);

% S.Valence.n = n; S.Valence.l = l; S.Valence.f_up = f_up; S.Valence.f_dw = f_dw;

% Only l = 0,1,2,3 possible
index0 = l==0; % selecting all l = 0
n0 = n(index0); % all 'n' for l = 0
f0_tot = f_tot(index0); % corresponding total occupation
f0_up = f_up(index0); % corresponding up-spin occupation
f0_dw = f_dw(index0); % corresponding dw-spin occupation

index1 = l==1; % selecting all l = 1
n1 = n(index1); % all 'n' for l = 1
f1_tot = f_tot(index1); % corresponding total occupation
f1_up = f_up(index1); % corresponding up-spin occupation
f1_dw = f_dw(index1); % corresponding dw-spin occupation

index2 = l==2; % selecting all l = 2
n2 = n(index2); % all 'n' for l = 2
f2_tot = f_tot(index2); % corresponding total occupation
f2_up = f_up(index2); % corresponding up-spin occupation
f2_dw = f_dw(index2); % corresponding dw-spin occupation

index3 = l==3; % selecting all l = 3 
n3 = n(index3); % all 'n' for l = 3
f3_tot = f_tot(index3); % corresponding total occupation
f3_up = f_up(index3); % corresponding up-spin occupation
f3_dw = f_dw(index3); % corresponding dw-spin occupation

if S.usefock < 2
    fprintf("===========================================================\n")
    valence = [n;l;f_tot];
    T = table(valence,'RowNames',{'n','l','Occupation'},'VariableNames',...
        {'Valence'});
    disp(T)
    fprintf("===========================================================\n")
end

S.n0 = n0; S.n1 = n1; S.n2 = n2; S.n3 = n3; S.f_up = f_up; S.f_dw = f_dw;
S.f0_up = f0_up; S.f1_up = f1_up; S.f2_up = f2_up; S.f3_up = f3_up;
S.f0_dw = f0_dw; S.f1_dw = f1_dw; S.f2_dw = f2_dw; S.f3_dw = f3_dw;
S.f0_tot = f0_tot; S.f1_tot = f1_tot; S.f2_tot = f2_tot; S.f3_tot = f3_tot;
%--------------------------------------------------------------------------
% First SCF loop with guess density

fprintf("\nStarting SCF loop with guess-density.\n\n")
fprintf("======================SCF Iteration #%d======================\n",1)

if S.xc == 4 % first step in mGGA is GGA
    fprintf("This step is GGA_PBE.\n")
end

% Poisson Solve
t_poiss = tic;
phi = poissonSolve(S, input_struct);
t_poiss = toc(t_poiss); 
fprintf("Time for Poisson-Solve:\t\t %0.10f s\n", t_poiss)

% XC
t_XC = tic;
S.countPotential = -1;
[S, Vxc] = xcPotentialAtom(S);
t_XC = toc(t_XC);
fprintf("Time for calculating V_xc:\t %0.10f s\n", t_XC)

orb_count = 1;
orbitals = zeros(Nd-1,2*length(n));
orbital_l = zeros(length(n),2);
eigen_values = zeros(length(n),2);
rho_up = 0;
rho_dw = 0;
e_occ_up = [];
e_occ_dw = [];
t_eigS = tic;

% Calculate rho according to occupation
for i = min(l):max(l)
    
    if i == 0
        prin_quantum_possible = n0;
        o_tot = f0_tot;
        o_up = f0_up;
        o_dw = f0_dw;
    elseif i == 1
        prin_quantum_possible = n1;
        o_tot = f1_tot;
        o_up = f1_up;
        o_dw = f1_dw;
    elseif i == 2
        prin_quantum_possible = n2;
        o_tot = f2_tot;
        o_up = f2_up;
        o_dw = f2_dw;
    else
        prin_quantum_possible = n3;
        o_tot = f3_tot;
        o_up = f3_up;
        o_dw = f3_dw;
    end
    
    len = length(prin_quantum_possible);
    [S, eigvec, eigval] = eigenSolve(phi, Vxc, i, S);
    
    for j = 1:len
        % Normalizing eigen-vectors (orbitals) : Spin-up
        func = [0;eigvec(:,j);0];
        integral = w*((func.^2)./int_scale);
        eigvec(:,j) = (1/sqrt(integral))*eigvec(:,j);
        
        % Normalizing eigen-vectors (orbitals) : Spin-dw
        func = [0;eigvec(:,j+(Nd-1));0];
        integral = w*((func.^2)./int_scale);
        eigvec(:,j+(Nd-1)) = (1/sqrt(integral))*eigvec(:,j+(Nd-1));
        
        % Store Occupation
        e_occ_up = [e_occ_up; o_up(j)];
        e_occ_dw = [e_occ_dw; o_dw(j)];
                
        % Calculate rho
        rho_up = rho_up + o_up(j)*(eigvec(:,j).^2)./(r(2:Nd).^2);
        rho_dw = rho_dw + o_dw(j)*(eigvec(:,j+(Nd-1)).^2)./(r(2:Nd).^2);
        
        orbitals(:,orb_count) = eigvec(:,j);
        orbital_l(orb_count,1) = i;
        orbitals(:,orb_count+length(n)) = eigvec(:,j+(Nd-1));
        orbital_l(orb_count,2) = i;
        
        eigen_values(orb_count,1) = eigval(j,1);
        eigen_values(orb_count,2) = eigval(j,2);
        orb_count = orb_count+1;
    end
end
t_eigS = toc(t_eigS); S.tEigS = S.tEigS+t_eigS;
fprintf("Time for Eigen-Solve:\t\t %0.10f s\n\n", t_eigS)
S.orbitals.matrix = orbitals;
S.orbital_l.matrix = orbital_l;
S.occ.matrix = [e_occ_up e_occ_dw];
rho_up = 0.25*rho_up/pi;
rho_dw = 0.25*rho_dw/pi;
rho = (rho_up + rho_dw);

iter = 1;
if ~S.spinFlag
    rho_prev = S.rho(:,1);  
else
    rho_prev = S.rho(:,2:3);
end

S.rho(:,1) = rho;
S.rho(:,2) = rho_up;
S.rho(:,3) = rho_dw;

% Initialize the mixing history vectors
S.X = zeros((Nd-1)*S.nspin,S.MixingHistory);
S.F = zeros((Nd-1)*S.nspin,S.MixingHistory);
S.mixing_hist_fkm1 = zeros(Nd-1,1);

if ~S.spinFlag
    S.mixing_hist_xkm1 = rho;
else
    S.mixing_hist_xkm1 = reshape(S.rho(:,2:3),[],1);
end

%--------------------------------------------------------------------------
while iter < S.MAXIT_SCF
    % Simple Density Mixing
    
    % Error
    if ~S.spinFlag
        error = norm(S.rho(:,1)-rho_prev)/norm(S.rho(:,1));
    else
        integrand = [0; (S.r(2:S.Nd).^2).*(S.mag);0];
        S.netM = 4*pi*S.w*(integrand./S.int_scale);
        fprintf('Net magnetization:\t\t ');
        fprintf('%.6f\n', S.netM);
        rho_updw = S.rho(:,2:3);
        error = norm(reshape(rho_updw - rho_prev,[],1))/norm(rho_updw(:));
    end
    
    fprintf('Relative Error in Density:\t %0.14e\n\n',error)
    
    % Periodic Pulay Mixing
    if ~S.spinFlag
        [S, S.rho(:,1)] = Periodic_Pulay(S, S.rho(:,1), rho_prev, iter, input_struct);
        rho_prev = S.rho(:,1);
    else
        [S, rho_updw(:)] = Periodic_Pulay(S, rho_updw(:), rho_prev(:), iter, input_struct);
        S.rho(:,2:3) = rho_updw;
        S.rho(:,1) = S.rho(:,2) + S.rho(:,3);
        S.mag = S.rho(:,2) - S.rho(:,3);
        rho_prev = S.rho(:,2:3);
    end
    
       
    if error > scf_tol
        fprintf("======================SCF Iteration #%d======================\n",iter+1)
        iter = iter+1;
        S.negrhoFlag = 0;
    elseif error < scf_tol && iter > 1 && S.negrhoFlag == 0
        fprintf("Density has already converged to %0.4e.\n",scf_tol)
        fprintf("Total number of SCF \t\t: %d \n",iter)
        break
    elseif error < scf_tol && iter == 1 
        fprintf("======================SCF Iteration #%d======================\n",iter+1)
        iter = iter+1;
    elseif error < scf_tol && S.negrhoFlag == 1
        fprintf("======================SCF Iteration #%d======================\n",iter+1)
        iter = iter+1;
        S.neghrhoFlag = 0;
    end
        
    % Poisson and Vxc Solve
    t_poiss = tic;
    phi = poissonSolve(S, input_struct);
    t_poiss = toc(t_poiss);
    fprintf("Time for Poisson-Solve:\t\t %0.10f s\n", t_poiss)
    
    t_XC = tic;
    [S, Vxc] = xcPotentialAtom(S);
    t_XC = toc(t_XC);
    fprintf("Time for calculating V_xc:\t %0.10f s\n", t_XC)
    
    % Re-initialize
    orb_count = 1;
    orbitals = zeros(Nd-1,2*length(n));
    eigen_values = zeros(length(n),2);
    rho_up = 0;
    rho_dw = 0;
    total_eig_time = 0;
    t_eigS = tic;
    % Calculate rho according to occupation
    for i = min(l):max(l)
        
        total_eig_time = total_eig_time+t_eigS;
        if i == 0
            prin_quantum_possible = n0;
            o_tot = f0_tot;
            o_up = f0_up;
            o_dw = f0_dw;
        elseif i == 1
            prin_quantum_possible = n1;
            o_tot = f1_tot;
            o_up = f1_up;
            o_dw = f1_dw;
        elseif i == 2
            prin_quantum_possible = n2;
            o_tot = f2_tot;
            o_up = f2_up;
            o_dw = f2_dw;
        else
            prin_quantum_possible = n3;
            o_tot = f3_tot;
            o_up = f3_up;
            o_dw = f3_dw;
        end
        
        len = length(prin_quantum_possible);
        [S, eigvec, eigval] = eigenSolve(phi, Vxc, i, S);
        
        for j = 1:len
            % Normalizing eigen-vectors (orbitals) : Spin-up
            func = [0;eigvec(:,j);0];
            integral = w*((func.^2)./int_scale);
            eigvec(:,j) = (1/sqrt(integral))*eigvec(:,j);
            
            % Normalizing eigen-vectors (orbitals) : Spin-dw
            func = [0;eigvec(:,j+(Nd-1));0];
            integral = w*((func.^2)./int_scale);
            eigvec(:,j+(Nd-1)) = (1/sqrt(integral))*eigvec(:,j+(Nd-1));
            
            % Calculate rho
            rho_up = rho_up + o_up(j)*(eigvec(:,j).^2)./(r(2:Nd).^2);
            rho_dw = rho_dw + o_dw(j)*(eigvec(:,j+(Nd-1)).^2)./(r(2:Nd).^2);
            
            orbitals(:,orb_count) = eigvec(:,j);
            orbitals(:,orb_count+length(n)) = eigvec(:,j+(Nd-1));
            
            eigen_values(orb_count,1) = eigval(j,1);
            eigen_values(orb_count,2) = eigval(j,2);
            orb_count = orb_count+1;
        end
        
    end
    t_eigS = toc(t_eigS);S.tEigS = S.tEigS+t_eigS;
    fprintf("Time for Eigen-Solve:\t\t %0.10f s\n\n", t_eigS)
    S.orbitals.matrix = orbitals;
    rho_up = 0.25*rho_up/pi;
    rho_dw = 0.25*rho_dw/pi;
    
    S.rho(:,1) = rho_up + rho_dw;
    S.rho(:,2) = rho_up;
    S.rho(:,3) = rho_dw;

end
S.EigVal = eigen_values;
S.phi = phi;
S.Vxc = Vxc;
if iter == S.MAXIT_SCF
    fprintf("<strong>The SCF cycle did not converge!</strong>\n")
%     fprintf("============================================================\n")
%     fprintf("************************************************************\n")
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Desnity Matrix
function D = densityMatrix(S)
%--------------------------------------------------------------------------
Nd = S.Nd;
w = S.w;
int_scale = S.int_scale;

orbitals = S.orbitals.matrix;
orbital_l = S.orbital_l.matrix(:);
occ = S.occ.matrix;
occ_up = occ(:,1); occ_dw = occ(:,2);
%--------------------------------------------------------------------------

lmin = min(orbital_l); lmax = max(orbital_l);
D = zeros(Nd+1,Nd+1);


for i = lmin : lmax
    if i == 0
        jstop = length(S.n0);
        jstart = 1;
    elseif i == 1
        jstop = length(S.n0)+length(S.n1);
        jstart = length(S.n0)+1;
    elseif i == 2
        jstop = length(S.n0)+length(S.n1)+length(S.n2);
        jstart = length(S.n0)+length(S.n1)+1;
    elseif i == 3
        jstop = length(S.n0)+length(S.n1)+length(S.n2)+length(S.n3);
        jstart = length(S.n0)+length(S.n1)+length(S.n2)+1;
    end
    
    for j = jstart : jstop
        orbital = orbitals(:,j);
        orbital = [0; orbital; 0];
        
        wt_orbital = (w'.*orbital)./int_scale;
        V = occ_up(j)*(orbital*wt_orbital');
        D = D + V;
    end
end

spin_up_last_col = length(S.n0)+length(S.n1)+length(S.n2)...
    +length(S.n3);
for i = lmin : lmax
    if i == 0
        jstop = spin_up_last_col+length(S.n0);
        jstart = spin_up_last_col+1;
    elseif i == 1
        jstop = spin_up_last_col+length(S.n0)+length(S.n1);
        jstart = length(S.n0)+spin_up_last_col+1;
    elseif i == 2
        jstop = spin_up_last_col+length(S.n0)+length(S.n1)+length(S.n2);
        jstart = length(S.n0)+length(S.n1)+spin_up_last_col+1;
    elseif i == 3
        jstop = spin_up_last_col+length(S.n0)+length(S.n1)+length(S.n2)+length(S.n3);
        jstart = length(S.n0)+length(S.n1)+length(S.n2)+spin_up_last_col+1;
    end
    
    for j = jstart : jstop
        orbital = orbitals(:,j);
        orbital = [0; orbital; 0];
        
        wt_orbital = (w'.*orbital)./int_scale;
        V = occ_dw(j-spin_up_last_col)*(orbital*wt_orbital');
        D = D + V;
    end
end

D = D(2:Nd,2:Nd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Store Exx operator
function S = exx_operator(S)
orbital_l = S.orbital_l.matrix(:);
lmin = min(orbital_l); lmax = max(orbital_l);

for l = lmin : lmax
    spin = 0.5;
    V = evaluateExxPotential(S,l,spin);
    S.l_channel(l+1).Vexx_up = V;
end
if S.spinFlag
    for l = lmin:lmax
        spin = -0.5;
        V = evaluateExxPotential(S,l,spin);
        S.l_channel(l+1).Vexx_dw = V;
    end
end
end