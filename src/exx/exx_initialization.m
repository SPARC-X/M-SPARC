function S = exx_initialization(S)

if S.MAXIT_FOCK < 0
    S.MAXIT_FOCK = 20;
end
if S.MINIT_FOCK < 0
    S.MINIT_FOCK = 2;
end
if S.FOCK_TOL < 0
    S.FOCK_TOL = 0.2* S.SCF_tol;
end
if S.SCF_tol_init < 0
    S.SCF_tol_init = max(10*S.FOCK_TOL,1e-3);
end

if S.xc == 40
    S.hyb_mixing = 1.0;
elseif S.xc == 41
    S.hyb_mixing = 0.25;
elseif S.xc == 427
    S.hyb_mixing = 0.25;
    if S.hyb_range_fock < 0
        S.hyb_range_fock = 0.1587;          % VASP
    end
    if S.hyb_range_pbe < 0
        S.hyb_range_pbe = 0.1587;           % VASP
    end
    % S.hyb_range_fock = 0.15/sqrt(2);  % ABINIT
    % S.hyb_range_pbe = 0.15*2^(1/3);   % ABINIT
    % S.hyb_range_fock = 0.106;         % QE
    % S.hyb_range_pbe = 0.106;          % QE
end
if strcmp(S.ExxMethod, 'FOURIER_SPACE') || strcmp(S.ExxMethod, 'fourier_space')
    S.exxmethod = 0;
elseif strcmp(S.ExxMethod, 'REAL_SPACE') || strcmp(S.ExxMethod, 'real_space')
    S.exxmethod = 1;
    if S.BC ~= 1
        error('Real space solver is only available in molecule simulation.\n');
    end
elseif strcmp(S.ExxMethod, '')
    fprintf("Default setting: Solving Exact Exchange in Fourier Space.\n");
    S.exxmethod = 0;
    S.ExxMethod = 'FOURIER_SPACE';
else
    error('Please provide correct method for solving Exact Exchange, fourier_space or real_space.\n');
end

if strcmp(S.ExxDivMethod, 'SPHERICAL') || strcmp(S.ExxDivMethod, 'spherical')
    S.exxdivmethod = 0;
elseif strcmp(S.ExxDivMethod, 'AUXILIARY') || strcmp(S.ExxDivMethod, 'auxiliary')
    S.exxdivmethod = 1;
elseif strcmp(S.ExxDivMethod, 'ERFC') || strcmp(S.ExxDivMethod, 'erfc')
    S.exxdivmethod = 2;
    if S.xc ~= 427
        error('Error: ERFC method could only be used with HSE functional.\n');
    end
elseif strcmp(S.ExxDivMethod, '')
    if S.BC <= 2
        if S.xc == 427
            fprintf("Default setting: Using simple erfc for singularity issue for bulk or molecule with HSE.\n");
            S.exxdivmethod = 2;
            S.ExxDivMethod = 'ERFC';
        else
            fprintf("Default setting: Using spherical truncation for singularity issue for bulk or molecule with hybrid functional.\n");
            S.exxdivmethod = 0;
            S.ExxDivMethod = 'SPHERICAL';
        end
    else
        fprintf("Default setting: Using auxiliary function method for singularity issue for slab and wire with hybrid functional.\n");
        S.exxdivmethod = 1;
        S.ExxDivMethod = 'AUXILIARY';
    end
else
    error('Please provide correct method for singularity in Exact Exchange, spherical, auxiliary or erfc.\n');
end

if S.isgamma == 0
    S.exxmethod = 0;    % only fourier space method is allowed
end
% ensure all occupied states are used. 
if S.EXXACEVal_state < 0
    S.EXXACEVal_state = 3;
end
if sum(S.exx_downsampling == [1,1,1]) ~= 3
    if sum(S.exx_downsampling < 0) > 0
        error('Please provide non-negative integer for EXX_DOWNSAMPLING.\n');
    end
    for i = 1:3
        if S.exx_downsampling(i) ~= 0
            if mod(S.nkpt(i),S.exx_downsampling(i))
                error("Number of kpoints should be divisible by EXX_DOWNSAMPLING.\n");
            end
        end
    end
end
S = const_for_FFT(S);
S = kshift_phasefactor(S);
end


function S = const_for_FFT(S)
N1 = S.Nx;
N2 = S.Ny;
N3 = S.Nz;

L1 = N1 * S.dx;
L2 = N2 * S.dy;
L3 = N3 * S.dz;

tnkpthf = S.tnkpthf;
tnkpt = S.tnkpt;
TOL = 1e-8;

gmet = S.grad_T * S.grad_T' ;
ucvol = L1*L2*L3*S.Jacb*S.tnkpthf;
R_c = (3*ucvol/(4*pi))^(1/3);
    
if S.cell_typ < 3
    sumx = 0;
    sumy = 0; 
    sumz = 0;
end

S.kpthf_ind = zeros(tnkpthf,2);
% 1 -> k, 0 -> -k
% find the index of kptgridhf in reduced kptgrid
% TODO: use ismembertol
for i = 1:tnkpthf
    for j = 1:tnkpt
        if (abs(S.kptgridhf(i,1) - S.kptgrid(j,1)) < TOL && abs(S.kptgridhf(i,2) - S.kptgrid(j,2)) < TOL && abs(S.kptgridhf(i,3) - S.kptgrid(j,3)) < TOL)
            S.kpthf_ind(i,1) = j;
            S.kpthf_ind(i,2) = 1;
            break;
        end
    end
    if S.kpthf_ind(i,2)
        continue;
    end
    for j = 1:tnkpt
        if (abs(S.kptgridhf(i,1) + S.kptgrid(j,1) - sumx) < TOL && abs(S.kptgridhf(i,2) + S.kptgrid(j,2) - sumy) < TOL && abs(S.kptgridhf(i,3) + S.kptgrid(j,3) - sumz) < TOL)
            S.kpthf_ind(i,1) = j;
            break;
        end
    end
end


% find unique kpoint shift
shift = zeros(tnkpthf*tnkpt,3);
count = 1;
for k_index = 1:tnkpt
    for q_index = 1:tnkpthf
        k = S.kptgrid(k_index,:);
        q = S.kptgridhf(q_index,:);
        shift(count,:) = k - q;
        count = count + 1;
    end
end
S.k_shift = uniquetol(shift,TOL,'ByRows',true);
S.num_shift = size(S.k_shift,1);
% always put zero shift in the end of list 
zero_ind = find(ismembertol(S.k_shift,[0,0,0],1e-8,'ByRows',true))+0;
S.k_shift([zero_ind,S.num_shift],:) = S.k_shift([S.num_shift,zero_ind],:);

S.const_by_alpha = zeros(S.num_shift,N1,N2,N3);

if S.exxdivmethod == 0
    if S.Calc_stress
        S.const_stress = zeros(S.num_shift,N1,N2,N3);
        S.const_stress_2 = zeros(S.num_shift,N1,N2,N3);
    end
    % spherical truncation method by Spencer 
    for ind = 1:S.num_shift
        count = 1;
        Gpkmq2 = zeros(N1,N2,N3);
        for k3 = [1:floor(N3/2)+1, floor(-N3/2)+2:0]
            for k2 = [1:floor(N2/2)+1, floor(-N2/2)+2:0]
                for k1 = [1:floor(N1/2)+1, floor(-N1/2)+2:0]
                    G = [(k1-1)*2*pi/L1, (k2-1)*2*pi/L2, (k3-1)*2*pi/L3];
                    Gpkmq = G + S.k_shift(ind,:);
                    Gpkmq2(count) = Gpkmq * (gmet * Gpkmq');
                    count = count + 1;
                end
            end
        end
        iszero = Gpkmq2 < 1e-4;
        Gpkmq2(iszero) = 1;
        const = 1 - cos(R_c*sqrt(Gpkmq2));
        const(iszero) = R_c^2/2;
        S.const_by_alpha(ind,:,:,:) = 4*pi*const./Gpkmq2;
        
        if S.Calc_stress
            x = R_c*sqrt(Gpkmq2);
            const = 1 - cos(x) - x/2.*sin(x);
            const(iszero) = R_c^4/24;
            S.const_stress(ind,:,:,:) = 4*pi*const./(Gpkmq2.^2);
            
            % 1/3 factor copied from ABINIT.
            const = 0.5*x.*sin(x);
            const(iszero) = R_c^2/2;
            S.const_stress_2(ind,:,:,:) = 4*pi*const./Gpkmq2/3;
        end
    end
    
elseif S.exxdivmethod == 1
    if S.Calc_stress
        S.const_stress = zeros(S.num_shift,N1,N2,N3);
    elseif S.Calc_pres
        S.const_press = zeros(S.num_shift,N1,N2,N3);
    end
    % auxiliary function method by Gygi
    aux = exx_divergence(S);
    for ind = 1:S.num_shift
        count = 1;
        Gpkmq2 = zeros(N1,N2,N3);
        for k3 = [1:floor(N3/2)+1, floor(-N3/2)+2:0]
            for k2 = [1:floor(N2/2)+1, floor(-N2/2)+2:0]
                for k1 = [1:floor(N1/2)+1, floor(-N1/2)+2:0]
                    G = [(k1-1)*2*pi/L1, (k2-1)*2*pi/L2, (k3-1)*2*pi/L3];
                    Gpkmq = G + S.k_shift(ind,:);
                    Gpkmq2(count) = Gpkmq * (gmet * Gpkmq');
                    count = count + 1;
                end
            end
        end
        iszero = Gpkmq2 < 1e-4;
        Gpkmq2(iszero) = 1;
        if S.hyb_range_fock > 0
            const = 1 - exp(-0.25/S.hyb_range_fock^2*Gpkmq2);
            const(iszero) = aux + 0.25/S.hyb_range_fock^2;
        else
            const = ones(N1,N2,N3);
            const(iszero) = aux;
        end

        S.const_by_alpha(ind,:,:,:) = 4*pi*const./Gpkmq2;
        
        if S.Calc_stress
            if S.hyb_range_fock > 0
                x = -0.25/S.hyb_range_fock^2*Gpkmq2;
                const = 1 - exp(x).*(1-x);
                const(iszero) = (0.25/S.hyb_range_fock^2)^2;
                S.const_stress(ind,:,:,:) = 4*pi*const./(Gpkmq2.^2);
            else
                const = ones(N1,N2,N3);
                const(iszero) = 0;
                S.const_stress(ind,:,:,:) = 4*pi*const./(Gpkmq2.^2);
            end
            
            % for consistency of stress formula
            S.const_stress(ind,:,:,:) = S.const_stress(ind,:,:,:)/4;
            
        elseif S.Calc_pres
            if S.hyb_range_fock > 0
                x = -0.25/S.hyb_range_fock^2*Gpkmq2;
                const = 1 - exp(x).*(1-x);
                const(iszero) = 0;
                S.const_press(ind,:,:,:) = 4*pi*const./Gpkmq2;
            else
                const = ones(N1,N2,N3);
                const(iszero) = 0;
                S.const_press(ind,:,:,:) = 4*pi*const./Gpkmq2;
            end
            
            % for consistency of pressure formula
            S.const_press(ind,:,:,:) = S.const_press(ind,:,:,:)/4;
        end
    end
    
elseif S.exxdivmethod == 2
    if S.Calc_stress
        S.const_stress = zeros(S.num_shift,N1,N2,N3);
    elseif S.Calc_pres
        S.const_press = zeros(S.num_shift,N1,N2,N3);
    end
    % Simple method by ERFC
    for ind = 1:S.num_shift
        count = 1;
        Gpkmq2 = zeros(N1,N2,N3);
        for k3 = [1:floor(N3/2)+1, floor(-N3/2)+2:0]
            for k2 = [1:floor(N2/2)+1, floor(-N2/2)+2:0]
                for k1 = [1:floor(N1/2)+1, floor(-N1/2)+2:0]
                    G = [(k1-1)*2*pi/L1, (k2-1)*2*pi/L2, (k3-1)*2*pi/L3];
                    Gpkmq = G + S.k_shift(ind,:);
                    Gpkmq2(count) = Gpkmq * (gmet * Gpkmq');
                    count = count + 1;
                end
            end
        end
        iszero = Gpkmq2 < 1e-4;
        Gpkmq2(iszero) = 1;
        const = 1 - exp(-0.25/S.hyb_range_fock^2*Gpkmq2);
        const(iszero) = 0.25/S.hyb_range_fock^2;
        
        S.const_by_alpha(ind,:,:,:) = 4*pi*const./Gpkmq2;
        
        if S.Calc_stress
            x = -0.25/S.hyb_range_fock^2*Gpkmq2;
            const = 1 - exp(x).*(1-x);
            const(iszero) = 0;
            S.const_stress(ind,:,:,:) = 4*pi*const./(Gpkmq2.^2);
        elseif S.Calc_pres
            x = -0.25/S.hyb_range_fock^2*Gpkmq2;
            const = 1 - exp(x).*(1-x);
            const(iszero) = 0;
            S.const_press(ind,:,:,:) = 4*pi*const./Gpkmq2;
        end
    end
end

end


function S = kshift_phasefactor(S)
if S.num_shift == 1
    return;
end
r = zeros(S.N,3);
count = 1;
for k = 0:S.Nz-1
    for j = 0:S.Ny-1
        for i = 0:S.Nx-1
            r(count,1) = i/S.Nx*S.L1;
            r(count,2) = j/S.Ny*S.L2;
            r(count,3) = k/S.Nz*S.L3;
            count = count + 1;
        end
    end
end

neg_phase = zeros(S.N,S.num_shift-1);
pos_phase = zeros(S.N,S.num_shift-1);
for i = 1:S.num_shift-1
    k_shift = S.k_shift(i,:);
    neg_phase(:,i) = exp(-1i*r*k_shift');
    pos_phase(:,i) = exp(1i*r*k_shift');
end
S.neg_phase = neg_phase;
S.pos_phase = pos_phase;
end




function c = exx_divergence(S)
N1 = S.Nx;
N2 = S.Ny;
N3 = S.Nz;

L1 = N1 * S.dx;
L2 = N2 * S.dy;
L3 = N3 * S.dz;
gmet = S.grad_T * S.grad_T';

% alpha is related to ecut
ecut = ecut_estimate(S.dx,S.dy,S.dz);
alpha = 10/(ecut*2);
fprintf(' Ecut estimation is %.2f Ha (%.2f Ry) and alpha in auxiliary function is %.6f\n',ecut, ecut*2, alpha);

sumfGq = 0;
% Use QE's choice of q vector
% Actually all q vectors by Monkhorst-pack grid also works

MPG_typ = @(nkpt) (0:nkpt-1); % MP grid points for finite group order

kptgrid_x = (1/S.nkpt(1)) * MPG_typ(S.nkpthf(1));
kptgrid_y = (1/S.nkpt(2)) * MPG_typ(S.nkpthf(2));
kptgrid_z = (1/S.nkpt(3)) * MPG_typ(S.nkpthf(3));

[kptgrid_X, kptgrid_Y, kptgrid_Z] = ndgrid(kptgrid_x,kptgrid_y,kptgrid_z);
kptgrid = [reshape(kptgrid_X,[],1),reshape(kptgrid_Y,[],1),reshape(kptgrid_Z,[],1)];
    
for q = 1:S.tnkpthf
%     xq = S.kptgridhf(khf,:)./[2*pi/L1,2*pi/L2,2*pi/L3];
    xq = kptgrid(q,:);
    
    for k3 = [1:floor(N3/2)+1, floor(-N3/2)+2:0]
        for k2 = [1:floor(N2/2)+1, floor(-N2/2)+2:0]
            for k1 = [1:floor(N1/2)+1, floor(-N1/2)+2:0]
                G = [(k1-1), (k2-1), (k3-1)];

                Gpq = (G+xq).*[2*pi/L1,2*pi/L2,2*pi/L3];
                Gpq2 = Gpq * (gmet * Gpq');
                if Gpq2 > 1e-8
                    if S.hyb_range_fock > 0
                        sumfGq = sumfGq + exp(-alpha*Gpq2)/Gpq2*(1-exp(-Gpq2/4/S.hyb_range_fock^2));
                    else
                        sumfGq = sumfGq + exp(-alpha*Gpq2)/Gpq2;
                    end
                end
            end
        end
    end
end

V = L1*L2*L3*S.Jacb;
Nk = S.tnkpthf;

if S.hyb_range_fock > 0
    sumfGq = sumfGq + 1/4/S.hyb_range_fock^2;
    
    nqq = 100000;
    dq = 5/sqrt(alpha)/nqq;
    aa = 0;
    for iq = 0:nqq
        q_ = dq*(iq+0.5);
        qq = q_*q_;
        aa = aa - exp( -alpha * qq) * exp(-qq/4/S.hyb_range_fock^2)*dq;
    end
    aa = aa*2/pi + 1/sqrt(pi*alpha);
    scaled_intf = V*Nk/(4*pi)*aa;
    c = scaled_intf - sumfGq;
else
    
    scaled_intf = V*Nk/(4*pi*sqrt(pi*alpha));
    c = scaled_intf + alpha - sumfGq;
end

fprintf(' The constant for zero G is %f\n',c);
% exx_div = -c*2*4*pi in qe
end


function ecut = ecut_estimate(hx,hy,hz)
dx2_inv = 1/(hx * hx);
dy2_inv = 1/(hy * hy);
dz2_inv = 1/(hz * hz);
h_eff = sqrt(3.0 / (dx2_inv + dy2_inv + dz2_inv));

ecut = (pi/h_eff)^2/2;                      % ecut_rho
ecut = ecut / 4;                            % ecut_psi
ecut = ecut * 0.9;                          % by experience
end


function S = const_for_FFT_fd(S)
N1 = S.Nx;
N2 = S.Ny;
N3 = S.Nz;

dx = S.dx;
dy = S.dy;
dz = S.dz;

L1 = N1 * S.dx;
L2 = N2 * S.dy;
L3 = N3 * S.dz;

w1 = S.w1;
w2 = S.w2;

dx2 = dx*dx; dy2 = dy*dy; dz2 = dz*dz;

T11 = S.lapc_T(1,1);T22 = S.lapc_T(2,2);T33 = S.lapc_T(3,3);
T12 = S.lapc_T(2,1);T23 = S.lapc_T(3,2);T13 = S.lapc_T(3,1);

w2_x = w2 * T11 / dx2;
w2_y = w2 * T22 / dy2;
w2_z = w2 * T33 / dz2;
w2_diag = w2_x(1) + w2_y(1) +w2_z(1);

w1_x = w1 / dx;
w1_y = w1 / dy;
w1_z = w1 / dz;

FDn = S.FDn;

tnkpthf = S.tnkpthf;
tnkpt = S.tnkpt;
TOL = 1e-8;

gmet = S.grad_T * S.grad_T' ;
ucvol = L1*L2*L3*S.Jacb*S.tnkpthf;
R_c = (3*ucvol/(4*pi))^(1/3);
    
if S.cell_typ < 3
    sumx = 0;
    sumy = 0; 
    sumz = 0;
end

S.kpthf_ind = zeros(tnkpthf,2);
% 1 -> k, 0 -> -k
% find the index of kptgridhf in reduced kptgrid
% TODO: use ismembertol
for i = 1:tnkpthf
    for j = 1:tnkpt
        if (abs(S.kptgridhf(i,1) - S.kptgrid(j,1)) < TOL && abs(S.kptgridhf(i,2) - S.kptgrid(j,2)) < TOL && abs(S.kptgridhf(i,3) - S.kptgrid(j,3)) < TOL)
            S.kpthf_ind(i,1) = j;
            S.kpthf_ind(i,2) = 1;
            break;
        end
    end
    if S.kpthf_ind(i,2)
        continue;
    end
    for j = 1:tnkpt
        if (abs(S.kptgridhf(i,1) + S.kptgrid(j,1) - sumx) < TOL && abs(S.kptgridhf(i,2) + S.kptgrid(j,2) - sumy) < TOL && abs(S.kptgridhf(i,3) + S.kptgrid(j,3) - sumz) < TOL)
            S.kpthf_ind(i,1) = j;
            break;
        end
    end
end


% find unique kpoint shift
shift = zeros(tnkpthf*tnkpt,3);
count = 1;
for k_index = 1:tnkpt
    for q_index = 1:tnkpthf
        k = S.kptgrid(k_index,:);
        q = S.kptgridhf(q_index,:);
        shift(count,:) = k - q;
        count = count + 1;
    end
end
S.k_shift = uniquetol(shift,TOL,'ByRows',true);
S.num_shift = size(S.k_shift,1);
% always put zero shift in the end of list 
zero_ind = find(ismembertol(S.k_shift,[0,0,0],1e-8,'ByRows',true))+0;
S.k_shift([zero_ind,S.num_shift],:) = S.k_shift([S.num_shift,zero_ind],:);

S.const_by_alpha = zeros(S.num_shift,N1,N2,N3);

if S.exxdivmethod == 0
    if S.Calc_stress
        S.const_stress = zeros(S.num_shift,N1,N2,N3);
        S.const_stress_2 = zeros(S.num_shift,N1,N2,N3);
    end
elseif S.exxdivmethod == 1
    aux = exx_divergence(S);
    if S.Calc_stress
        S.const_stress = zeros(S.num_shift,N1,N2,N3);
    elseif S.Calc_pres
        S.const_press = zeros(S.num_shift,N1,N2,N3);
    end
elseif S.exxdivmethod == 2
    if S.Calc_stress
        S.const_stress = zeros(S.num_shift,N1,N2,N3);
    elseif S.Calc_pres
        S.const_press = zeros(S.num_shift,N1,N2,N3);
    end
end

for ind = 1:S.num_shift
    kpt1 = S.k_shift(ind,1);
    kpt2 = S.k_shift(ind,2);
    kpt3 = S.k_shift(ind,3);
    count = 1;
    Gpkmq2 = zeros(N1,N2,N3);
    for k3 = [1:floor(N3/2)+1, floor(-N3/2)+2:0]
        for k2 = [1:floor(N2/2)+1, floor(-N2/2)+2:0]
            for k1 = [1:floor(N1/2)+1, floor(-N1/2)+2:0]
                % finite difference G2
		        Gpkmq2(count) = -w2_diag;
                for p = 1:FDn
                    Gpkmq2(count) = Gpkmq2(count) - 2 * ...
                        (  cos(2*pi*(k1-1)*p/N1+kpt1*dx*p)*w2_x(p+1) ...
                         + cos(2*pi*(k2-1)*p/N2+kpt2*dy*p)*w2_y(p+1) ...
                         + cos(2*pi*(k3-1)*p/N3+kpt3*dz*p)*w2_z(p+1));
                end
                if S.cell_typ == 2
                    for p = 1:FDn
                        for q = 1:FDn
                            Gpkmq2(count) = Gpkmq2(count) + 8 *...
                                T12 * w1_x(p+1) * w1_y(q+1) * sin(2*pi*(k1-1)*p/N1+kpt1*dx*p) * sin(2*pi*(k2-1)*q/N2+kpt2*dy*q);
                            Gpkmq2(count) = Gpkmq2(count) + 8 *...
                                T13 * w1_x(p+1) * w1_z(q+1) * sin(2*pi*(k1-1)*p/N1+kpt1*dx*p) * sin(2*pi*(k3-1)*q/N3+kpt3*dz*q);
                            Gpkmq2(count) = Gpkmq2(count) + 8 *...
                                T23 * w1_y(p+1) * w1_z(q+1) * sin(2*pi*(k2-1)*p/N2+kpt2*dy*p) * sin(2*pi*(k3-1)*q/N3+kpt3*dz*q);
                        end
                    end
                end
                count = count + 1;
            end
        end
    end

    if S.exxdivmethod == 0
        iszero = Gpkmq2 < 1e-4;
        Gpkmq2(iszero) = 1;
        const = 1 - cos(R_c*sqrt(Gpkmq2));
        const(iszero) = R_c^2/2;
        S.const_by_alpha(ind,:,:,:) = 4*pi*const./Gpkmq2;
        
        if S.Calc_stress
            x = R_c*sqrt(Gpkmq2);
            const = 1 - cos(x) - x/2.*sin(x);
            const(iszero) = R_c^4/24;
            S.const_stress(ind,:,:,:) = 4*pi*const./(Gpkmq2.^2);
            
            % 1/3 factor copied from ABINIT.
            const = 0.5*x.*sin(x);
            const(iszero) = R_c^2/2;
            S.const_stress_2(ind,:,:,:) = 4*pi*const./Gpkmq2/3;
        end
    elseif S.exxdivmethod == 1
        iszero = Gpkmq2 < 1e-4;
        Gpkmq2(iszero) = 1;
        if S.hyb_range_fock > 0
            const = 1 - exp(-0.25/S.hyb_range_fock^2*Gpkmq2);
            const(iszero) = aux + 0.25/S.hyb_range_fock^2;
        else
            const = ones(N1,N2,N3);
            const(iszero) = aux;
        end
        S.const_by_alpha(ind,:,:,:) = 4*pi*const./Gpkmq2;
        
        if S.Calc_stress
            if S.hyb_range_fock > 0
                x = -0.25/S.hyb_range_fock^2*Gpkmq2;
                const = 1 - exp(x).*(1-x);
                const(iszero) = (0.25/S.hyb_range_fock^2)^2;
                S.const_stress(ind,:,:,:) = 4*pi*const./(Gpkmq2.^2);
            else
                const = ones(N1,N2,N3);
                const(iszero) = 0;
                S.const_stress(ind,:,:,:) = 4*pi*const./(Gpkmq2.^2);
            end
            
            % for consistency of stress formula
            S.const_stress(ind,:,:,:) = S.const_stress(ind,:,:,:)/4;
            
        elseif S.Calc_pres
            if S.hyb_range_fock > 0
                x = -0.25/S.hyb_range_fock^2*Gpkmq2;
                const = 1 - exp(x).*(1-x);
                const(iszero) = 0;
                S.const_press(ind,:,:,:) = 4*pi*const./Gpkmq2;
            else
                const = ones(N1,N2,N3);
                const(iszero) = 0;
                S.const_press(ind,:,:,:) = 4*pi*const./Gpkmq2;
            end
            
            % for consistency of pressure formula
            S.const_press(ind,:,:,:) = S.const_press(ind,:,:,:)/4;
        end
    elseif S.exxdivmethod == 2
        iszero = Gpkmq2 < 1e-4;
        Gpkmq2(iszero) = 1;
        const = 1 - exp(-0.25/S.hyb_range_fock^2*Gpkmq2);
        const(iszero) = 0.25/S.hyb_range_fock^2;
        S.const_by_alpha(ind,:,:,:) = 4*pi*const./Gpkmq2;
        if S.Calc_stress
            x = -0.25/S.hyb_range_fock^2*Gpkmq2;
            const = 1 - exp(x).*(1-x);
            const(iszero) = 0;
            S.const_stress(ind,:,:,:) = 4*pi*const./(Gpkmq2.^2);
        elseif S.Calc_pres
            x = -0.25/S.hyb_range_fock^2*Gpkmq2;
            const = 1 - exp(x).*(1-x);
            const(iszero) = 0;
            S.const_press(ind,:,:,:) = 4*pi*const./Gpkmq2;
        end
    end
end
end