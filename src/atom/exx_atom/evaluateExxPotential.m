function VexxL = evaluateExxPotential(S, l, spin)
% @brief    Calculates the exact exchange potential for each 'l' channel.
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
r = S.r;
int_scale = S.int_scale;

orbitals = S.SCF_orbitals.matrix;
orbital_l = S.orbital_l.matrix(:);
occ = S.occ.matrix;
occ_up = occ(:,1); occ_dw = occ(:,2);
%--------------------------------------------------------------------------

lmin = min(orbital_l); lmax = max(orbital_l);

VexxL_up = zeros(Nd+1,Nd+1);
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
        V = zeros(Nd+1,Nd+1);
        
        for k = abs(l - i) :2: abs(l + i)
            wt_orbital = (w'.*orbital)./int_scale;
            term = (orbital*wt_orbital')*(Wigner3j(l,i,k)^2);
            [row, col] = size(term);
            
            for rows = 1:row
                for cols = 1:col
                    ratio = ((min(r(rows),r(cols)))^k) / ((max(r(rows),r(cols)))^(k+1));
                    term(rows,cols) = ratio*term(rows,cols);
                end
            end
            
            V = V + term;
        end
        
        V = occ_up(j)*V;
        VexxL_up = VexxL_up + V;
    end
end

VexxL_dw = zeros(Nd+1,Nd+1);
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
    
    for j = jstart : jstop
        orbital = orbitals(:,j);
        orbital = [0; orbital; 0];
        V = zeros(Nd+1,Nd+1);
        
        for k = abs(l - i) :2: abs(l + i)
            wt_orbital = (w'.*orbital)./int_scale;
            term = (orbital*wt_orbital')*(Wigner3j(l,i,k)^2);
            [row, col] = size(term);
            
            for rows = 1:row
                for cols = 1:col
                    ratio = ((min(r(rows),r(cols)))^k) / ((max(r(rows),r(cols)))^(k+1));
                    term(rows,cols) = ratio*term(rows,cols);
                end
            end
            
            V = V + term;
        end
        
        V = occ_dw(j-spin_up_last_col)*V;
        VexxL_dw = VexxL_dw + V;
    end
end

if ~S.spinFlag
    VexxL = VexxL_up + VexxL_dw;
    VexxL = -0.5*VexxL;
else
    if spin == 0.5
        VexxL = -VexxL_up;
    else
        VexxL = -VexxL_dw;
    end
end

VexxL = VexxL(2:Nd,2:Nd);
end

