function [DX_x,DX_y,DX_z] = dpseudopot(X,II,JJ,KK,XX,YY,ZZ,S)
% @ brief      Function to perform gradient of a 3d matrix.
% @ authors
%         Abhiraj Sharma <asharma424@gatech.edu>
%         Qimen Xu <qimenxu@gatech.edu>
%         Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%====================================================================== 
    if S.cell_typ <= 2
        DX_x = zeros(size(X)); DX_y = zeros(size(X)); DX_z = zeros(size(X));
        for p = 1:S.FDn
            DX_x(II,JJ,KK) = DX_x(II,JJ,KK) + S.w1(p+1)/S.dx*(X(II+p,JJ,KK)-X(II-p,JJ,KK));
            DX_y(II,JJ,KK) = DX_y(II,JJ,KK) + S.w1(p+1)/S.dy*(X(II,JJ+p,KK)-X(II,JJ-p,KK));
            DX_z(II,JJ,KK) = DX_z(II,JJ,KK) + S.w1(p+1)/S.dz*(X(II,JJ,KK+p)-X(II,JJ,KK-p));
        end
    elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
        [DX_x,DX_y,DX_z] = dpseudopot_cychel(X,II,JJ,KK,XX,YY,ZZ,S);
    end
end