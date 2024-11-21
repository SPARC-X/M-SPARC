function [S, Vxc] = xcPotentialAtom(S)
% @brief    Calculates the XC potential and energy for the atom case.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Boqin Zhang <bzhang376@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param S      Struct containing relevant fields
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================

% rho: total density
% sigma: square of gradient of rho
% exc: the energy per unit particle
% vxc: first partial derivative of the energy per unit volume in terms of the density
% v2xc: first partial derivative of the energy per unit volume in terms of sigma

if S.spinFlag == 0
    rho = S.rho(:,1);
    if S.NLCC_flag
        rho = rho+S.rho_Tilde;
    end
    rho(rho < S.xc_rhotol) = S.xc_rhotol;
    Nd = S.Nd;
    r = S.r(2:Nd);
    if S.isGradient
        D = S.Gradient.matrix;
        D = D(2:Nd,2:Nd);
        rho_r = r.*rho;       
        drho = (D*rho_r - rho)./r;
        
        sigma = drho.*drho;
        sigma(sigma < S.xc_rhotol) = S.xc_rhotol;
        S.sigma = sigma;
    end
    
    % metaGGA
    if S.ixc(3) == 1 
        if S.countPotential == -1 % compute potential before SCF by GGA_PBE, there is no psi yet
            S.ixc = [2 3 1 0];
            S.xc_option = [1 1];
        else % compute kinetic energy density tau
            [S, S.tau] = kineticDensity(S, rho); % it is \tau=-1/2\nabla^2\phi = 1/2(\nabla\phi)\cdot(\nabla\phi)
        end
    end
    
    % iexch
    switch S.ixc(1)
        case 1
            [ex, vx] = slater(rho);
            v2x = zeros(size(rho));
        case 2
            [ex,vx,v2x] = pbex(rho,sigma,S.xc_option(1));
        case 3
            [ex,vx,v2x] = rPW86x(rho,sigma);
        case 4
            [ex,vx,v2x,v3x] = scanx(rho,sigma,S.tau);
        case 5
            [ex,vx,v2x,v3x] = rscanx(rho,sigma,S.tau);
        case 6
            [ex,vx,v2x,v3x] = r2scanx(rho,sigma,S.tau);
        otherwise
            ex = zeros(size(rho));
            vx = zeros(size(rho));
            v2x = zeros(size(rho));
    end
    
    % icorr
    switch S.ixc(2)
        case 1
            [ec, vc] = pz(rho);
            v2c = zeros(size(rho));
        case 2
            [ec, vc] = pw(rho);
            v2c = zeros(size(rho));
        case 3
            [ec,vc,v2c] = pbec(rho,sigma,S.xc_option(2));
        case 4
            [ec,vc,v2c,v3c] = scanc(rho,sigma,S.tau);
        case 5
            [ec,vc,v2c,v3c] = rscanc(rho,sigma,S.tau); 
        case 6
            [ec,vc,v2c,v3c] = r2scanc(rho,sigma,S.tau);
        otherwise
            ec = zeros(size(rho));
            vc = zeros(size(rho));
            v2c = zeros(size(rho));
    end
    
    % hybrid
    if S.usefock > 1
        if S.xc == 41
            ex = ex*(1 - S.hyb_mixing);
            vx = vx*(1 - S.hyb_mixing);
            v2x = v2x*(1 - S.hyb_mixing);
        end
    end
    
    % imeta
    if S.ixc(3) == 1
        if S.countPotential == -1 % compute potential before SCF by GGA_PBE, there is no psi yet
            if strcmp(S.XC, 'SCAN')
                S.ixc = [4 4 1 0]; % restore the labels
            elseif strcmp(S.XC, 'RSCAN')
                S.ixc = [5 5 1 0]; % restore the labels
            else % R2SCAN
                S.ixc = [6 6 1 0]; % restore the labels
            end
        else
            S.VxcScan3 = v3x + v3c;
        end
        S.countPotential = S.countPotential + 1;
    end
    
    
    exc = ex + ec;
    vxc = vx + vc;
    v2xc = v2x + v2c;
    
    if S.isGradient
        Vxc = vxc - D*(v2xc.*drho) - (2./r).*v2xc.*drho;
    else
        Vxc = vxc;
    end
else
    rho = S.rho;
    if S.NLCC_flag
        rho(:,2) = rho(:,2)+S.rho_Tilde * 0.5;
        rho(:,3) = rho(:,3)+S.rho_Tilde * 0.5;
    end
    rho(rho < S.xc_rhotol) = S.xc_rhotol;
    rho(:,1) = rho(:,2) + rho(:,3);
    Nd = S.Nd;
    r = S.r(2:Nd);
    
    if S.isGradient
        D = S.Gradient.matrix;
        D = D(2:Nd,2:Nd);
        rho_tot_r = r.*rho(:,1);
        rho_up_r = r.*rho(:,2);
        rho_dw_r = r.*rho(:,3);
        
        drho_tot= (D*rho_tot_r - rho(:,1))./r;
        drho_up = (D*rho_up_r - rho(:,2))./r;
        drho_dw = (D*rho_dw_r - rho(:,3))./r;
        drho = [drho_tot, drho_up, drho_dw];
        
        sigma_tot = drho_tot.*drho_tot;
        sigma_up = drho_up.*drho_up;
        sigma_dw = drho_dw.*drho_dw;
        sigma = [sigma_tot, sigma_up, sigma_dw];
        sigma(sigma < S.xc_rhotol) = S.xc_rhotol;
        S.sigma = sigma;
    end
    
    if S.ixc(3) == 1 % metaGGA
        if S.countPotential == -1 % compute potential before SCF by GSGA_PBE, there is no psi yet
            S.ixc = [2 3 1 0];
            S.xc_option = [1 1];
        else % compute kinetic energy density tau
            [S S.tau] = kineticDensity(S, rho); % it is \tau=-1/2\nabla^2\phi = 1/2(\nabla\phi)\cdot(\nabla\phi)
        end
    end
    
    % iexch
    switch S.ixc(1)
        case 1
            [ex, vx] = slater_spin(rho);
            v2x = zeros(Nd-2,2);
        case 2
            [ex, vx, v2x] = pbex_spin(rho,sigma,S.xc_option(1));
        case 3
            [ex,vx,v2x] = rPW86x_spin(rho,sigma);
        case 4
            [ex,vx,v2x,v3x] = scanx_spin(rho,sigma,S.tau);
        case 5
            [ex,vx,v2x,v3x] = rscanx_spin(rho,sigma,S.tau);
        case 6
            [ex,vx,v2x,v3x] = r2scanx_spin(rho,sigma,S.tau);
        otherwise
            ex = zeros(Nd-2,1);
            vx = zeros(Nd-2,2);
            v2x = zeros(Nd-2,2);
    end
    
    % icorr
    switch S.ixc(2)
        case 1
            [ec, vc] = pz_spin(rho);
            v2c = zeros(Nd-2,1);
        case 2
            [ec, vc] = pw_spin(rho);
            v2c = zeros(Nd-2,1);
        case 3
            [ec, vc, v2c] = pbec_spin(rho,sigma,S.xc_option(1));
        case 4
            [ec,vc,v2c,v3c] = scanc_spin(rho,sigma,S.tau);
        case 5
            [ec,vc,v2c,v3c] = rscanc_spin(rho,sigma,S.tau); 
        case 6
            [ec,vc,v2c,v3c] = r2scanc_spin(rho,sigma,S.tau);
        otherwise
            ec = zeros(Nd-2,1);
            vc = zeros(Nd-2,2);
            v2c = zeros(Nd-2,1);
    end
    
    % hybrid
    if S.usefock > 1
        if S.xc == 41
            ex = ex*(1 - S.hyb_mixing);
            vx = vx*(1 - S.hyb_mixing);
            v2x = v2x*(1 - S.hyb_mixing);
        end
    end
    
    exc = ex + ec;
    vxc = vx + vc;
    v2xc = [v2c v2x];
    
    % imeta
    if S.ixc(3) == 1
        if S.countPotential == -1 % compute potential before SCF by GGA_PBE, there is no psi yet
            if strcmp(S.XC, 'SCAN')
                S.ixc = [4 4 1 0]; % restore the labels
            elseif strcmp(S.XC, 'RSCAN')
                S.ixc = [5 5 1 0]; % restore the labels
            else % R2SCAN
                S.ixc = [6 6 1 0]; % restore the labels
            end
        else
            S.VxcScan3 = v3x + v3c;
        end
        S.countPotential = S.countPotential + 1;
    end
    
    
    if S.isGradient
        Vxc_temp = D*(v2xc.*drho) + (2./r).*v2xc.*drho;
        vxc(:,1) = vxc(:,1) - Vxc_temp(:,2) - Vxc_temp(:,1);
        vxc(:,2) = vxc(:,2) - Vxc_temp(:,3) - Vxc_temp(:,1);
    end
    Vxc = vxc;
end

S.exc = exc;

end

%--------------------------------------------------------------------------
function [S, tau] = kineticDensity(S, rho)
% Eq A2 of Sala, Fabiano and Constantin, Phys. Rev. B, 91, 035126 (2015)
% (2l + 1) factor disappears due to equal smearing
% Eq 17 of Lehtola,S.; J. Chem. Theory Comput. 2023, 19, 2502−2517.

tau = zeros(size(rho,1), size(rho,2)); 
Nd = S.Nd;
D = S.Gradient.matrix;
D = D(2:Nd,2:Nd);
r = S.r(2:Nd);
occ = S.occ.matrix(:)';
R_nl_Tilde = (S.orbitals.matrix);
R_nl = R_nl_Tilde./r;
orbital_l = S.orbital_l.matrix(:);
dR_nl_dr = ((D*R_nl_Tilde) - R_nl)./r;
dR_nl_dr_sq = dR_nl_dr.*dR_nl_dr;


for spinor = 1:S.nspinor
    if S.spinFlag == 0
        start = 1;
        stop = size(R_nl,2);
    else
        if spinor == 1
            start = 1;
            stop = size(R_nl,2)/2;
        else
            start = size(R_nl,2)/2 + 1;
            stop = size(R_nl,2);
        end
    end
    
    for i = start : stop
        l = orbital_l(i);
        term1 = 0.5*occ(i)*dR_nl_dr_sq(:,i);
        term2 = 0.5*occ(i)*(l*(l+1))*((R_nl(:,i).^2)./(r.^2));
        term = (term1 + term2);
        tau(:,spinor) = tau(:,spinor) + term;
    end
end

tau = (0.25/pi)*tau;

if size(tau,2) == 1
    tau = tau;
else
    tau = [sum(tau(:,1:2),2) tau(:,1:2)];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LDA XC - No Spin
function [ex,vx] = slater(rho)
% slater exchange
% @param rho  = total electron density

% parameter 
C2 = 0.73855876638202; % 3/4*(3/pi)^(1/3)
C3 = 0.9847450218427;  % (3/pi)^(1/3)

% computation
ex = - C2 * rho.^(1./3.);
vx = - C3 * rho.^(1./3.);
end

function [ex,vx] = slater_spin(rho)
% slater exchange, spin-polarized case
% @param rho  = [rho_tot rho_up rho_dw]

% parameters 
third = 1./3.;
threefourth_divpi = 3.0/4.0/pi;
sixpi2_1_3 = (6.0 * pi^2)^third;

% computationns
rho_tot = rho(:,1);
rho_updw = rho(:,2:3);
rho_updnm1_3 = rho_updw.^(-third);
rhom1_3 = rho_tot.^(-third);
rhotot_inv = rhom1_3.^3;
ex_lsd = -threefourth_divpi * sixpi2_1_3 * (rho_updnm1_3 .* rho_updnm1_3 .* rho_updw);

vx = (4/3) * ex_lsd;
exc = sum(ex_lsd .* rho_updw,2);
ex = exc .* rhotot_inv;
end

function [ec,vc] = pz(rho)
% pz correlation
% @param rho  = total electron density
% J.P. Perdew and A. Zunger, PRB 23, 5048 (1981).

% parameters
A = 0.0311;
B = -0.048;
C = 0.002;
D = -0.0116;
gamma1 = -0.1423;
beta1 = 1.0529;
beta2 = 0.3334;

% compuatation
ec = zeros(size(rho,1),1);
vc = zeros(size(rho,1),1);
rs = (0.75./(pi*rho)).^(1.0/3.0) ;
islt1 = (rs < 1.0);
lnrs = log(rs(islt1));
sqrtrs = sqrt(rs(~islt1));
ec(islt1) = A * lnrs + B + C * rs(islt1) .* lnrs + D * rs(islt1);
ox = 1.0 + beta1*sqrtrs + beta2*rs(~islt1);
ec(~islt1) = gamma1 ./ ox;
vc(islt1) = lnrs.*(A + (2.0/3.0)*C*rs(islt1)) + (B-(1.0/3.0)*A) + (1.0/3.0)*(2.0*D-C)* rs(islt1);
vc(~islt1) = ec(~islt1) .* (1 + (7.0/6.0)*beta1*sqrtrs + (4.0/3.0)*beta2*rs(~islt1)) ./ ox;
end

function [ec,vc] = pw(rho)
% pw correlation
% @param rho  = total electron density
% J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)

% parameters
A = 0.031091 ;
alpha1 = 0.21370 ;
beta1 = 7.5957 ;
beta2 = 3.5876 ;
beta3 = 1.6382 ;
beta4 = 0.49294 ;

% computation
rs = (0.75./(pi*rho)).^(1./3.);
rsm12 = rs.^(-0.5);
rs12 = rs.^0.5;
rs32 = rs.^1.5;
rs2 = rs.^2;

om = 2*A*(beta1*rs12 + beta2*rs + beta3*rs32 + beta4*rs2);
dom = A*(beta1*rsm12+ 2*beta2 + 3*beta3*rs12 + 2*2*beta4*rs);
olog = log(1 + 1./om);
t = -2*A*(1+alpha1*rs);
ec = t.*olog;
vc = ec - (rs/3.).*(-2*A*alpha1*olog - (t.*dom)./(om.*om+om) ) ;
end

function [ex,v1x,v2x] = rPW86x(rho,sigma)
a = 1.851;
b = 17.33;
c = 0.163;
s_prefactor = 6.18733545256027; % 2*(3\pi^2)^(1/3)
Ax = -0.738558766382022; % -3/4 * (3/pi)^(1/3)
four_thirds = 4.0/3.0;

grad_rho = sigma.^0.5;
s = grad_rho ./ (s_prefactor*rho.^four_thirds);
s_2 = s.*s;
s_3 = s_2.*s;
s_4 = s_3.*s;
s_5 = s_3.*s_2;
s_6 = s_5.*s;
fs = (1.0 + a*s_2 + b*s_4 + c*s_6).^(1.0/15.0);
ex = Ax * rho.^(1.0/3.0) .* fs; % \epsilon_x, not n\epsilon_x
df_ds = (1.0./(15.0*fs.^14.0)) .* (2.0*a*s + 4.0*b*s_3 + 6.0*c*s_5);
v1x = Ax*four_thirds * (rho.^(1.0/3.0) .*fs - grad_rho./(s_prefactor*rho).*df_ds);
v2x = Ax * df_ds./(s_prefactor*grad_rho);
end

function [Ec, Vc] = LDA_VWNc(rho)
for i = 1:length(rho)
    [ec,vc] = LDA_VWNc_scalar(rho(i));
    Vc(i,1) = vc;
    Ec(i,1) = ec;
end
end

function [ec, vc] = LDA_VWNc_scalar(rho)
% Used to get vxc for a scalar rho

if (abs(rho) < 1e-14)
    ex =0; ec = 0; vx = 0; vc = 0;
    return
end
% Parameters
y0 = -0.10498;
b = 3.72744;
c = 12.9352;
A = 0.0621814;

Q = sqrt(4*c - b^2);
rs = (3/(4*pi*rho))^(1/3);
y = sqrt(rs);

ec = A/2 * (log(y^2/get_Y(y, b, c)) + 2*b/Q * atan(Q/(2*y+b))  ...
   - b*y0/get_Y(y0, b, c) * ( ...
            log((y-y0)^2 / get_Y(y, b, c)) ...
            + 2*(b+2*y0) / Q * atan(Q/(2*y+b)) ...
          ) );
vc = ec - A/6 * (c*(y-y0)-b*y0*y)/((y-y0)*get_Y(y, b, c));


    function Y = get_Y(y,b,c)
        Y = y^2 + b*y + c;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LDA XC - Spin
function [ec,vc] = pz_spin(rho)
error("ERROR: pz_spin has not been implemented yet.");
end

function [ec,vc] = pw_spin(rho)
% pw correlation, spin-polarized case
% @param rho  = total electron density
% @param zeta  = zeta = (rho_up - rho_dw) / rho_tot

% parameters
third = 1./3.;
alpha_zeta2 = 1.0 - 1.0e-6; % ABINIT
alpha_zeta = 1.0 - 1.0e-6; % ABINIT
% XC.alpha_zeta2 = 1.0; XC.alpha_zeta = 1.0; %LIBXC
rsfac = 0.6203504908994000;
sq_rsfac = sqrt(rsfac);
sq_rsfac_inv = 1.0/sq_rsfac;
ec0_aa = 0.031091; ec1_aa = 0.015545; mac_aa = 0.016887; % ABINIT
%ec0_aa = 0.0310907; ec1_aa = 0.01554535; mac_aa = 0.0168869; % LIBXC
ec0_a1 = 0.21370;  ec1_a1 = 0.20548;  mac_a1 = 0.11125;
ec0_b1 = 7.5957;  ec1_b1 = 14.1189;  mac_b1 = 10.357;
ec0_b2 = 3.5876;   ec1_b2 = 6.1977;   mac_b2 = 3.6231;
ec0_b3 = 1.6382;   ec1_b3 = 3.3662;   mac_b3 = 0.88026;
ec0_b4 = 0.49294;  ec1_b4 = 0.62517;  mac_b4 = 0.49671;
factf_zeta = 1.0/(2.0^(4.0/3.0) - 2.0);
factfp_zeta = 4.0/3.0 * factf_zeta * alpha_zeta2;
fsec_inv = 1.0/1.709921;

% computation
rho_tot = rho(:,1);
rhom1_3 = rho_tot.^(-third);
rhotot_inv = rhom1_3.^3;
zeta = (rho(:,2) - rho(:,3)) .* rhotot_inv; % Check whether it is rho_up-rho_dn or rho_dn-rho_up
zetp = 1 + zeta * alpha_zeta;
zetm = 1 - zeta * alpha_zeta;
zetpm1_3 = zetp.^(-third);
zetmm1_3 = zetm.^(-third);

rhotmo6 = sqrt(rhom1_3);
rhoto6 = rho_tot .* rhom1_3 .* rhom1_3 .* rhotmo6;

% -----------------------------------------------------------------------------
% Then takes care of the LSD correlation part of the functional
rs = rsfac * rhom1_3;
sqr_rs = sq_rsfac * rhotmo6;
rsm1_2 = sq_rsfac_inv * rhoto6;

% Formulas A6-A8 of PW92LSD
ec0_q0 = -2.0 * ec0_aa * (1.0 + ec0_a1 * rs);
ec0_q1 = 2.0 * ec0_aa *(ec0_b1 * sqr_rs + ec0_b2 * rs + ec0_b3 * rs .* sqr_rs + ec0_b4 * rs .* rs);
ec0_q1p = ec0_aa * (ec0_b1 * rsm1_2 + 2.0 * ec0_b2 + 3.0 * ec0_b3 * sqr_rs + 4.0 * ec0_b4 * rs);
ec0_den = 1.0./(ec0_q1 .* ec0_q1 + ec0_q1);
ec0_log = -log(ec0_q1 .* ec0_q1 .* ec0_den);
ecrs0 = ec0_q0 .* ec0_log;
decrs0_drs = -2.0 * ec0_aa * ec0_a1 * ec0_log - ec0_q0 .* ec0_q1p .* ec0_den;

mac_q0 = -2.0 * mac_aa * (1.0 + mac_a1 * rs);
mac_q1 = 2.0 * mac_aa * (mac_b1 * sqr_rs + mac_b2 * rs + mac_b3 * rs .* sqr_rs + mac_b4 * rs .* rs);
mac_q1p = mac_aa * (mac_b1 * rsm1_2 + 2 * mac_b2 + 3 * mac_b3 * sqr_rs + 4 * mac_b4 * rs);
mac_den = 1.0./(mac_q1 .* mac_q1 + mac_q1);
mac_log = -log( mac_q1 .* mac_q1 .* mac_den );
macrs = mac_q0 .* mac_log;
dmacrs_drs = -2.0 * mac_aa * mac_a1 * mac_log - mac_q0 .* mac_q1p .* mac_den;

ec1_q0 = -2.0 * ec1_aa * (1.0 + ec1_a1 * rs);
ec1_q1 = 2.0 * ec1_aa * (ec1_b1 * sqr_rs + ec1_b2 * rs + ec1_b3 * rs .* sqr_rs + ec1_b4 * rs .* rs);
ec1_q1p = ec1_aa * (ec1_b1 * rsm1_2 + 2 * ec1_b2 + 3 * ec1_b3 * sqr_rs + 4 * ec1_b4 * rs);
ec1_den = 1.0./(ec1_q1 .* ec1_q1 + ec1_q1);
ec1_log = -log( ec1_q1 .* ec1_q1 .* ec1_den );
ecrs1 = ec1_q0 .* ec1_log;
decrs1_drs = -2.0 * ec1_aa * ec1_a1 * ec1_log - ec1_q0 .* ec1_q1p .* ec1_den;

% alpha_zeta is introduced in order to remove singularities for fully polarized systems.
zetp_1_3 = (1.0 + zeta * alpha_zeta) .* (zetpm1_3.^2);
zetm_1_3 = (1.0 - zeta * alpha_zeta) .* (zetmm1_3.^2);

f_zeta = ( (1.0 + zeta * alpha_zeta2) .* zetp_1_3 + (1.0 - zeta * alpha_zeta2) .* zetm_1_3 - 2.0 ) * factf_zeta;
fp_zeta = ( zetp_1_3 - zetm_1_3 ) * factfp_zeta;
zeta4 = zeta.^4;

gcrs = ecrs1 - ecrs0 + macrs * fsec_inv;
ecrs = ecrs0 + f_zeta .* (zeta4 .* gcrs - macrs * fsec_inv);
dgcrs_drs = decrs1_drs - decrs0_drs + dmacrs_drs * fsec_inv;
decrs_drs = decrs0_drs + f_zeta .* (zeta4 .* dgcrs_drs - dmacrs_drs * fsec_inv);
dfzeta4_dzeta = 4.0 * zeta.^3 .* f_zeta + fp_zeta .* zeta4;
decrs_dzeta = dfzeta4_dzeta .* gcrs - fp_zeta .* macrs * fsec_inv;
vxcadd = ecrs - rs * third .* decrs_drs - zeta .* decrs_dzeta;

ec = ecrs;
vc = zeros(size(rho,1),2);
vc(:,1) = vxcadd + decrs_dzeta;
vc(:,2) = vxcadd - decrs_dzeta;
end

function [ex,v1x,v2x] = rPW86x_spin(rho,sigma)
a = 1.851;
b = 17.33;
c = 0.163;
s_prefactor = 6.18733545256027; % 2*(3\pi^2)^(1/3)
Ax = -0.738558766382022; % -3/4 * (3/pi)^(1/3)
four_thirds = 4.0/3.0;

grad_rho = sigma(:,2:3).^0.5;
s = grad_rho ./ (2^(1/3) * s_prefactor*rho(:, 2:3).^four_thirds);
s_2 = s.*s;
s_3 = s_2.*s;
s_4 = s_3.*s;
s_5 = s_3.*s_2;
s_6 = s_5.*s;
fs = (1.0 + a*s_2 + b*s_4 + c*s_6).^(1.0/15.0);
ex = Ax*2^(1/3) * sum(rho(:, 2:3).^(1.0/3.0 + 1) .* fs, 2) ./ rho(:, 1); % \epsilon_x, not n\epsilon_x
df_ds = (1.0./(15.0*fs.^14.0)) .* (2.0*a*s + 4.0*b*s_3 + 6.0*c*s_5);
v1x = Ax*four_thirds * (2^(1/3)*rho(:, 2:3).^(1.0/3.0) .*fs - grad_rho./(s_prefactor*rho(:, 2:3)).*df_ds);
v2x = Ax * df_ds./(s_prefactor*grad_rho);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PBE XC - No Spin
function [ex,v1x,v2x] = pbex(rho,sigma,iflag)
% pbe exchange:
% @param rho  = total electron density
% @param grho = |\nabla rho|^2
% @param iflag options
% iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
% iflag=2  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
% iflag=3  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
% iflag=4  Zhang-Yang Revised PBE: Y. Zhang and W. Yang., Phys. Rev. Lett. 80, 890 (1998)
assert(iflag == 1 || iflag == 2 || iflag == 3 || iflag == 4);

% parameters 
mu_ = [0.2195149727645171 10.0/81.0 0.2195149727645171  0.2195149727645171];
mu = mu_(iflag);
kappa_ = [0.804 0.804 0.804 1.245];
kappa = kappa_(iflag);
threefourth_divpi = 3.0/4.0/pi;
sixpi2_1_3 = (6.0 * pi^2)^(1.0/3.0);
sixpi2m1_3 = 1.0/sixpi2_1_3;
mu_divkappa = mu/kappa;

% computation
rho_updn = rho/2.0;
rho_updnm1_3 = rho_updn.^(-1.0/3.0);
rhomot = rho_updnm1_3;
ex_lsd = -threefourth_divpi * sixpi2_1_3 * (rhomot .* rhomot .* rho_updn);
rho_inv = rhomot .* rhomot .* rhomot;
coeffss = (1.0/4.0) * sixpi2m1_3 * sixpi2m1_3 * (rho_inv .* rho_inv .* rhomot .* rhomot);
ss = (sigma/4.0) .* coeffss;

if iflag == 1 || iflag == 2 || iflag == 4
    divss = 1.0./(1.0 + mu_divkappa * ss);
    dfxdss = mu * (divss .* divss);
elseif iflag == 3
    divss = exp(-mu_divkappa * ss);
    dfxdss = mu * divss;
end

fx = 1.0 + kappa * (1.0 - divss);
dssdn = (-8.0/3.0) * (ss .* rho_inv);
dfxdn = dfxdss .* dssdn;
dssdg = 2.0 * coeffss;
dfxdg = dfxdss .* dssdg;

ex = ex_lsd .* fx;
v1x = ex_lsd .* ((4.0/3.0) * fx + rho_updn .* dfxdn);
v2x = 0.5 * ex_lsd .* rho_updn .* dfxdg;
end

function [ec,v1c,v2c] = pbec(rho,sigma,iflag)
% pbe correlation 
% @param rho  = total electron density
% @param sigma = |\nabla rho|^2
% iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
% iflag=2  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
% iflag=3  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
assert(iflag == 1 || iflag == 2 || iflag == 3);

% parameter 
beta_ = [0.066725 0.046 0.066725];
beta = beta_(iflag);
rsfac = 0.6203504908994000; % (0.75/pi)^(1/3)
sq_rsfac = sqrt(rsfac);
sq_rsfac_inv = 1.0/sq_rsfac;
third = 1.0/3.0;
twom1_3 = 2.0^(-third);
ec0_aa = 0.031091; 
ec0_a1 = 0.21370;  
ec0_b1 = 7.5957;
ec0_b2 = 3.5876;   
ec0_b3 = 1.6382;   
ec0_b4 = 0.49294;  
gamma = (1.0 - log(2.0)) /pi^2;
gamma_inv = 1/gamma;
coeff_tt = 1.0/(4.0 * 4.0 / pi * (3.0 * pi^2)^third);

% computation
rho_updn = rho/2.0;
rho_updnm1_3 = rho_updn.^(-third);
rhom1_3 = twom1_3 * rho_updnm1_3;

rhotot_inv = rhom1_3 .* rhom1_3 .* rhom1_3;
rhotmo6 = sqrt(rhom1_3);
rhoto6 = rho .* rhom1_3 .* rhom1_3 .* rhotmo6;

rs = rsfac * rhom1_3;
sqr_rs = sq_rsfac * rhotmo6;
rsm1_2 = sq_rsfac_inv * rhoto6;

%        Formulas A6-A8 of PW92LSD
ec0_q0 = -2.0 * ec0_aa * (1.0 + ec0_a1 * rs);
ec0_q1 = 2.0 * ec0_aa *(ec0_b1 * sqr_rs + ec0_b2 * rs + ec0_b3 * rs .* sqr_rs + ec0_b4 * rs .* rs);
ec0_q1p = ec0_aa * (ec0_b1 * rsm1_2 + 2.0 * ec0_b2 + 3.0 * ec0_b3 * sqr_rs + 4.0 * ec0_b4 * rs);
ec0_den = 1.0./(ec0_q1 .* ec0_q1 + ec0_q1);
ec0_log = -log(ec0_q1 .* ec0_q1 .* ec0_den);
ecrs0 = ec0_q0 .* ec0_log;
decrs0_drs = -2.0 * ec0_aa * ec0_a1 * ec0_log - ec0_q0 .* ec0_q1p .* ec0_den;

%        Add LSD correlation functional to GGA exchange functional
ec = ecrs0;
v1c = ecrs0 - (rs/3.0) .* decrs0_drs;

%        -----------------------------------------------------------------------------
%        Eventually add the GGA correlation part of the PBE functional
%        Note : the computation of the potential in the spin-unpolarized
%        case could be optimized much further. Other optimizations are left to do.

%        From ec to bb
bb = ecrs0 * gamma_inv;
dbb_drs = decrs0_drs * gamma_inv;

%        From bb to cc
exp_pbe = exp(-bb);
cc = 1.0./(exp_pbe - 1.0);
dcc_dbb = cc .* cc .* exp_pbe;
dcc_drs = dcc_dbb .* dbb_drs;

%        From cc to aa
coeff_aa = beta * gamma_inv;
aa = coeff_aa * cc;
daa_drs = coeff_aa * dcc_drs;

%        Introduce tt : do not assume that the spin-dependent gradients are collinear
dtt_dg = 2.0 * rhotot_inv .* rhotot_inv .* rhom1_3 * coeff_tt;
%        Note that tt is (the t variable of PBE divided by phi) squared
tt = 0.5 * sigma .* dtt_dg;

%        Get xx from aa and tt
xx = aa .* tt;
dxx_drs = daa_drs .* tt;
dxx_dtt = aa;

%        From xx to pade
pade_den = 1.0./(1.0 + xx .* (1.0 + xx));
pade = (1.0 + xx) .* pade_den;
dpade_dxx = -xx .* (2.0 + xx) .* (pade_den.^2);
dpade_drs = dpade_dxx .* dxx_drs;
dpade_dtt = dpade_dxx .* dxx_dtt;

%        From pade to qq
qq = tt .* pade;
dqq_drs = tt .* dpade_drs;
dqq_dtt = pade + tt .* dpade_dtt;

%        From qq to rr
arg_rr = 1.0 + beta * gamma_inv * qq;
div_rr = 1.0./arg_rr;
rr = gamma * log(arg_rr);
drr_dqq = beta * div_rr;
drr_drs = drr_dqq .* dqq_drs;
drr_dtt = drr_dqq .* dqq_dtt;

%        The GGA correlation energy is added
ec = ec + rr;

%        From hh to the derivative of the energy wrt the density
drhohh_drho = rr - third * rs .* drr_drs - (7.0/3.0) * tt .* drr_dtt; %- zeta * dhh_dzeta 
v1c = v1c + drhohh_drho;

%        From hh to the derivative of the energy wrt to the gradient of the
%        density, divided by the gradient of the density
%        (The v3.3 definition includes the division by the norm of the gradient)

v2c = rho .* dtt_dg .* drr_dtt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PBE XC - Spin
function [ex,v1x,v2x] = pbex_spin(rho,sigma,iflag)
% pbe exchange, spin-polarized case
% @param rho  = total electron density
% @param grho = |\nabla rho|^2
% @param iflag options
% iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
% iflag=2  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
% iflag=3  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
% iflag=4  Zhang-Yang Revised PBE: Y. Zhang and W. Yang., Phys. Rev. Lett. 80, 890 (1998)
assert(iflag == 1 || iflag == 2 || iflag == 3 || iflag == 4);

% parameters
third = 1./3.;
mu_ = [0.2195149727645171 10.0/81.0 0.2195149727645171 0.2195149727645171];
mu = mu_(iflag);
kappa_ = [0.804 0.804 0.804 1.245];
kappa = kappa_(iflag);
threefourth_divpi = 3.0/4.0/pi;
sixpi2_1_3 = (6.0 * pi^2)^(1.0/3.0);
sixpi2m1_3 = 1.0/sixpi2_1_3;
mu_divkappa = mu/kappa;

% computation
rho_tot = rho(:,1);
rho_updw = rho(:,2:3);
rho_updnm1_3 = rho_updw.^(-third);
rhom1_3 = rho_tot.^(-third);
rhotot_inv = rhom1_3.^3;

% -----------------------------------------------------------------------
% First take care of the exchange part of the functional

rhomot = rho_updnm1_3;
ex_lsd = -threefourth_divpi * sixpi2_1_3 * (rhomot .* rhomot .* rho_updw);

rho_inv = rhomot .* rhomot .* rhomot;
coeffss = (1.0/4.0) * sixpi2m1_3 * sixpi2m1_3 * (rho_inv .* rho_inv .* rhomot .* rhomot);
ss = sigma(:,2:3) .* coeffss;

if iflag == 1 || iflag == 2 || iflag == 4
    divss = 1.0./(1.0 + mu_divkappa * ss);
    dfxdss = mu * (divss .* divss);
elseif iflag == 3
    divss = exp(-mu_divkappa * ss);
    dfxdss = mu * divss;
end

fx = 1.0 + kappa * (1.0 - divss);
ex_gga = ex_lsd .* fx;
dssdn = (-8.0/3.0) * (ss .* rho_inv);
dfxdn = dfxdss .* dssdn;
v1x = ex_lsd .* ((4.0/3.0) * fx + rho_updw .* dfxdn);

dssdg = 2.0 * coeffss;
dfxdg = dfxdss .* dssdg; 
v2x = ex_lsd .* rho_updw .* dfxdg;
ex = sum(ex_gga .* rho_updw,2).* rhotot_inv;
end


function [ec,v1c,v2c] = pbec_spin(rho,sigma,iflag)
% pbe correlation, spin-polarized case
% @param rho  = total electron density
% @param sigma = |\nabla rho|^2
% iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
% iflag=2  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
% iflag=3  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
assert(iflag == 1 || iflag == 2 || iflag == 3);

% parameters
beta_ = [0.066725 0.046 0.066725];
beta = beta_(iflag);
rsfac = 0.6203504908994000; % (0.75/pi)^(1/3)
sq_rsfac = sqrt(rsfac);
sq_rsfac_inv = 1.0/sq_rsfac;
third = 1.0/3.0;
ec0_aa = 0.031091; ec1_aa = 0.015545; mac_aa = 0.016887; % ABINIT
%ec0_aa = 0.0310907; ec1_aa = 0.01554535; mac_aa = 0.0168869; % LIBXC
ec0_a1 = 0.21370;  ec1_a1 = 0.20548;  mac_a1 = 0.11125;
ec0_b1 = 7.5957;  ec1_b1 = 14.1189;  mac_b1 = 10.357;
ec0_b2 = 3.5876;   ec1_b2 = 6.1977;   mac_b2 = 3.6231;
ec0_b3 = 1.6382;   ec1_b3 = 3.3662;   mac_b3 = 0.88026;
ec0_b4 = 0.49294;  ec1_b4 = 0.62517;  mac_b4 = 0.49671;
gamma = (1.0 - log(2.0)) /pi^2;
gamma_inv = 1/gamma;
coeff_tt = 1.0/(4.0 * 4.0 / pi * (3.0 * pi^2)^third);
alpha_zeta = 1.0 - 1.0e-6; % ABINIT
% alpha_zeta = 1.0; %LIBXC
% alpha_zeta is introduced in order to remove singularities for fully polarized systems.
factf_zeta = 1.0/(2.0^(4.0/3.0) - 2.0);
factfp_zeta = 4.0/3.0 * factf_zeta * alpha_zeta;
fsec_inv = 1.0/1.709921;

% computation
rho_tot = rho(:,1);
rhom1_3 = rho_tot.^(-third);
rhotot_inv = rhom1_3.^3;
zeta = (rho(:,2) - rho(:,3)) .* rhotot_inv;
zetp = 1.0 + zeta * alpha_zeta;
zetm = 1.0 - zeta * alpha_zeta;
zetpm1_3 = zetp.^(-third);       
zetmm1_3 = zetm.^(-third);

rhotmo6 = sqrt(rhom1_3);
rhoto6 = rho_tot .* rhom1_3 .* rhom1_3 .* rhotmo6;

% -----------------------------------------------------------------------------
% Then takes care of the LSD correlation part of the functional

rs = rsfac * rhom1_3;
sqr_rs = sq_rsfac * rhotmo6;
rsm1_2 = sq_rsfac_inv * rhoto6;

%        Formulas A6-A8 of PW92LSD
ec0_q0 = -2.0 * ec0_aa * (1.0 + ec0_a1 * rs);
ec0_q1 = 2.0 * ec0_aa *(ec0_b1 * sqr_rs + ec0_b2 * rs + ec0_b3 * rs .* sqr_rs + ec0_b4 * rs .* rs);
ec0_q1p = ec0_aa * (ec0_b1 * rsm1_2 + 2.0 * ec0_b2 + 3.0 * ec0_b3 * sqr_rs + 4.0 * ec0_b4 * rs);
ec0_den = 1.0./(ec0_q1 .* ec0_q1 + ec0_q1);
ec0_log = -log(ec0_q1 .* ec0_q1 .* ec0_den);
ecrs0 = ec0_q0 .* ec0_log;
decrs0_drs = -2.0 * ec0_aa * ec0_a1 * ec0_log - ec0_q0 .* ec0_q1p .* ec0_den;

mac_q0 = -2.0 * mac_aa * (1.0 + mac_a1 * rs);
mac_q1 = 2.0 * mac_aa * (mac_b1 * sqr_rs + mac_b2 * rs + mac_b3 * rs .* sqr_rs + mac_b4 * rs .* rs);
mac_q1p = mac_aa * (mac_b1 * rsm1_2 + 2 * mac_b2 + 3 * mac_b3 * sqr_rs + 4 * mac_b4 * rs);
mac_den = 1.0./(mac_q1 .* mac_q1 + mac_q1);
mac_log = -log( mac_q1 .* mac_q1 .* mac_den );
macrs = mac_q0 .* mac_log;
dmacrs_drs = -2.0 * mac_aa * mac_a1 * mac_log - mac_q0 .* mac_q1p .* mac_den;

ec1_q0 = -2.0 * ec1_aa * (1.0 + ec1_a1 * rs);
ec1_q1 = 2.0 * ec1_aa * (ec1_b1 * sqr_rs + ec1_b2 * rs + ec1_b3 * rs .* sqr_rs + ec1_b4 * rs .* rs);
ec1_q1p = ec1_aa * (ec1_b1 * rsm1_2 + 2 * ec1_b2 + 3 * ec1_b3 * sqr_rs + 4 * ec1_b4 * rs);
ec1_den = 1.0./(ec1_q1 .* ec1_q1 + ec1_q1);
ec1_log = -log( ec1_q1 .* ec1_q1 .* ec1_den );
ecrs1 = ec1_q0 .* ec1_log;
decrs1_drs = -2.0 * ec1_aa * ec1_a1 * ec1_log - ec1_q0 .* ec1_q1p .* ec1_den;

zetp_1_3 = (1.0 + zeta * alpha_zeta) .* (zetpm1_3.^2);
zetm_1_3 = (1.0 - zeta * alpha_zeta) .* (zetmm1_3.^2);

f_zeta = ( (1.0 + zeta * alpha_zeta) .* zetp_1_3 + (1.0 - zeta * alpha_zeta) .* zetm_1_3 - 2.0 ) * factf_zeta;
fp_zeta = ( zetp_1_3 - zetm_1_3 ) * factfp_zeta;
zeta4 = zeta.^4;

gcrs = ecrs1 - ecrs0 + macrs * fsec_inv;
ecrs = ecrs0 + f_zeta .* (zeta4 .* gcrs - macrs * fsec_inv);

dgcrs_drs = decrs1_drs - decrs0_drs + dmacrs_drs * fsec_inv;
decrs_drs = decrs0_drs + f_zeta .* (zeta4 .* dgcrs_drs - dmacrs_drs * fsec_inv);
dfzeta4_dzeta = 4.0 * zeta.^3 .* f_zeta + fp_zeta .* zeta4;
decrs_dzeta = dfzeta4_dzeta .* gcrs - fp_zeta .* macrs * fsec_inv;

% Add LSD correlation functional to GGA exchange functional
ec = ecrs;
vxcadd = ecrs - rs * third .* decrs_drs - zeta .* decrs_dzeta;
v1c(:,1) = vxcadd + decrs_dzeta;
v1c(:,2) = vxcadd - decrs_dzeta;
% -----------------------------------------------------------------------------
% Eventually add the GGA correlation part of the PBE functional

% The definition of phi has been slightly changed, because
% the original PBE one gives divergent behaviour for fully polarized points

phi_zeta = ( zetpm1_3 .* (1.0 + zeta * alpha_zeta) + zetmm1_3 .* (1.0 - zeta * alpha_zeta)   ) * 0.5;
phip_zeta = (zetpm1_3 - zetmm1_3) * third * alpha_zeta;
phi_zeta_inv = 1.0./phi_zeta;
phi_logder = phip_zeta .* phi_zeta_inv;
phi3_zeta = phi_zeta .* phi_zeta .* phi_zeta;
gamphi3inv = gamma_inv * phi_zeta_inv .* phi_zeta_inv .* phi_zeta_inv;

%        From ec to bb
bb = ecrs .* gamphi3inv;
dbb_drs = decrs_drs .* gamphi3inv;
dbb_dzeta = gamphi3inv .* (decrs_dzeta - 3.0 * ecrs .* phi_logder);

% From bb to cc
exp_pbe = exp(-bb);
cc = 1.0./(exp_pbe - 1.0);
dcc_dbb = cc .* cc .* exp_pbe;
dcc_drs = dcc_dbb .* dbb_drs;
dcc_dzeta = dcc_dbb .* dbb_dzeta;

% From cc to aa
coeff_aa = beta * gamma_inv * phi_zeta_inv .* phi_zeta_inv;
aa = coeff_aa .* cc;
daa_drs = coeff_aa .* dcc_drs;
daa_dzeta = -2.0 * aa .* phi_logder + coeff_aa .* dcc_dzeta;

% Introduce tt : do not assume that the spin-dependent gradients are collinear
grrho2 = sigma(:,1);
dtt_dg = 2.0 * rhotot_inv .* rhotot_inv .* rhom1_3 * coeff_tt;
% Note that tt is (the t variable of PBE divided by phi) squared
tt = 0.5 * grrho2 .* dtt_dg;

% Get xx from aa and tt
xx = aa .* tt;
dxx_drs = daa_drs .* tt;
dxx_dzeta = daa_dzeta .* tt;
dxx_dtt = aa;

% From xx to pade
pade_den = 1.0./(1.0 + xx .* (1.0 + xx));
pade = (1.0 + xx) .* pade_den;
dpade_dxx = -xx .* (2.0 + xx) .* (pade_den.^2);
dpade_drs = dpade_dxx .* dxx_drs;
dpade_dtt = dpade_dxx .* dxx_dtt;
dpade_dzeta = dpade_dxx .* dxx_dzeta;

% From pade to qq
coeff_qq = tt .* phi_zeta_inv .* phi_zeta_inv;
qq = coeff_qq .* pade;
dqq_drs = coeff_qq .* dpade_drs;
dqq_dtt = pade .* phi_zeta_inv .* phi_zeta_inv + coeff_qq .* dpade_dtt;
dqq_dzeta = coeff_qq .* (dpade_dzeta - 2.0 * pade .* phi_logder);

% From qq to rr
arg_rr = 1.0 + beta * gamma_inv * qq;
div_rr = 1.0./arg_rr;
rr = gamma * log(arg_rr);
drr_dqq = beta * div_rr;
drr_drs = drr_dqq .* dqq_drs;
drr_dtt = drr_dqq .* dqq_dtt;
drr_dzeta = drr_dqq .* dqq_dzeta;

% From rr to hh
hh = phi3_zeta .* rr;
dhh_drs = phi3_zeta .* drr_drs;
dhh_dtt = phi3_zeta .* drr_dtt;
dhh_dzeta = phi3_zeta .* (drr_dzeta + 3.0 * rr .* phi_logder);
% The GGA correlation energy is added
ec = ec + hh;

% From hh to the derivative of the energy wrt the density
drhohh_drho = hh - third * rs .* dhh_drs - zeta .* dhh_dzeta - (7.0/3.0) * tt .* dhh_dtt;
v1c(:,1) = v1c(:,1) + drhohh_drho + dhh_dzeta;
v1c(:,2) = v1c(:,2) + drhohh_drho - dhh_dzeta;


% From hh to the derivative of the energy wrt to the gradient of the
% density, divided by the gradient of the density
% (The v3.3 definition includes the division by the norm of the gradient)

v2c = rho_tot .* dtt_dg .* dhh_dtt;           
end