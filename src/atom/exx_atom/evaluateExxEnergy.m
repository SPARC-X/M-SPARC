function S = evaluateExxEnergy(S)
% @brief    Calculates the exact exchange energy for the atom case.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param S      Struct containing relevant fields
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================
%--------------------------------------------------------------------------
Nd = S.Nd;
w = S.w;
int_scale = S.int_scale;

orbitals = S.SCF_orbitals.matrix;
orbital_l = S.orbital_l.matrix(:);
occ = S.occ.matrix;
occ_up = occ(:,1); occ_dw = occ(:,2);
lmin = min(orbital_l); lmax = max(orbital_l);

%--------------------------------------------------------------------------
Exx = 0;
% Exx_dens = zeros(size(S.rho));

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
    
    spin = 0.5;
    VexxL = evaluateExxPotential(S, i, spin);
    for j = jstart : jstop
        orbital = orbitals(:,j);
        wt_orbital = (w(2:Nd)'.*orbital)./int_scale(2:Nd);
        Exx = Exx + occ_up(j)*((wt_orbital')*(VexxL*orbital));
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
        jstop = spin_up_last_col+length(S.n0)+length(S.n1)...
            +length(S.n2);
        jstart = length(S.n0)+length(S.n1)+spin_up_last_col+1;
    elseif i == 3
        jstop = spin_up_last_col+length(S.n0)+length(S.n1)...
            +length(S.n2)+length(S.n3);
        jstart = length(S.n0)+length(S.n1)+length(S.n2)...
            +spin_up_last_col+1;
    end
    
    spin = -0.5;
    VexxL = evaluateExxPotential(S, i, spin);
    for j = jstart : jstop
        orbital = orbitals(:,j);
        wt_orbital = (w(2:Nd)'.*orbital)./int_scale(2:Nd);
        Exx = Exx + occ_dw(j-spin_up_last_col)*((wt_orbital')*(VexxL*orbital));
    end
end

Exx = 0.5*Exx;
S.Eex = S.hyb_mixing*Exx;
end

