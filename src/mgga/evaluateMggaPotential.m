function Hx = evaluateMggaPotential(S,X,Hx,kptvec,spin)
% @file    evaluateMggaPotential.m
% @brief   This file contains the functions applying Mgga potential to
%          Kohn-sham orbitals
% @authors Boqin Zhang <bzhang376@gatech.edu>
%          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>

if S.nspin == 1
    VxcScan3 = S.VxcScan3;
else % spin polarization mGSGA
    VxcScan3 = S.VxcScan3(:, spin); 
end
nCol = size(X, 2);
if S.cell_typ == 2 % unorthogonal cell
    lapc_T = [S.lapc_T(1,1), S.lapc_T(2,1), S.lapc_T(3,1);
        S.lapc_T(2,1), S.lapc_T(2,2), S.lapc_T(3,2);
        S.lapc_T(3,1), S.lapc_T(3,2), S.lapc_T(3,3)];
    v3grad1 = blochGradient(S,kptvec,1) *X; 
    v3grad2 = blochGradient(S,kptvec,2) *X; 
    v3grad3 = blochGradient(S,kptvec,3) *X; 

    v3gradpsiTheKpt = [v3grad1(:), v3grad2(:), v3grad3(:)];
    v3gradpsiTheKptMLapT = v3gradpsiTheKpt*lapc_T;
    v3gradpsiTheKptMLapT_1 = VxcScan3 .* reshape(v3gradpsiTheKptMLapT(:, 1), S.N, nCol);
    v3gradpsiTheKptMLapT_2 = VxcScan3 .* reshape(v3gradpsiTheKptMLapT(:, 2), S.N, nCol);
    v3gradpsiTheKptMLapT_3 = VxcScan3 .* reshape(v3gradpsiTheKptMLapT(:, 3), S.N, nCol);
    Hx = Hx - 0.5*(blochGradient(S,kptvec,1)*v3gradpsiTheKptMLapT_1...
                    + blochGradient(S,kptvec,2)*v3gradpsiTheKptMLapT_2...
                    + blochGradient(S,kptvec,3)*v3gradpsiTheKptMLapT_3);
else % orthogonal cell
    v3grad1 = VxcScan3 .* (blochGradient(S,kptvec,1) *X);
    v3grad2 = VxcScan3 .* (blochGradient(S,kptvec,2) *X);
    v3grad3 = VxcScan3 .* (blochGradient(S,kptvec,3) *X);
    Hx = Hx - 0.5*(blochGradient(S,kptvec,1)*v3grad1 ...
                    + blochGradient(S,kptvec,2)*v3grad2 ...
                    + blochGradient(S,kptvec,3)*v3grad3);
end

end