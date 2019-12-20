function b = pseudochargeDensity_atom(V,II,JJ,KK,xin,S)
% @ brief      Calculate lap * V = b (in this function)
% @ authors
%         Abhiraj Sharma <asharma424@gatech.edu>
%         Qimen Xu <qimenxu@gatech.edu>
%         Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param V           potential at FD nodes
% @param II          Node indices in X direction
% @param JJ          Node indices in Y direction
% @param KK          Node indices in Z direction
% @param xin         Reference x-coordinate
% @param b           laplacian of the potential
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%======================================================================     
	if S.cell_typ == 1
		b = zeros(size(V));
		dx2 = S.dx*S.dx;
		dy2 = S.dy*S.dy;
		dz2 = S.dz*S.dz;
		coeff =  S.w2(1) * (1/dx2 + 1/dy2 + 1/dz2);
		b(II,JJ,KK) = coeff * V(II,JJ,KK);
		for p = 1:S.FDn
			b(II,JJ,KK) = b(II,JJ,KK) + S.w2(p+1)/dx2 * (V(II+p,JJ,KK) + V(II-p,JJ,KK)) + ...
					S.w2(p+1)/dy2 * (V(II,JJ+p,KK) + V(II,JJ-p,KK)) + ...
					S.w2(p+1)/dz2 * (V(II,JJ,KK+p) + V(II,JJ,KK-p));
		end
	elseif S.cell_typ == 2
		b = zeros(size(V));
		dx2 = S.dx*S.dx;
		dy2 = S.dy*S.dy;
		dz2 = S.dz*S.dz;
		dxdy = S.dx*S.dy;
		dydz = S.dy*S.dz;
		dzdx = S.dz*S.dx;
		coeff =  S.w2(1) * (S.lapc_T(1,1)/dx2 + S.lapc_T(2,2)/dy2 + S.lapc_T(3,3)/dz2);
		b(II,JJ,KK) = coeff * V(II,JJ,KK);
		for p = 1:S.FDn
			b(II,JJ,KK) = b(II,JJ,KK) + S.w2(p+1)*S.lapc_T(1,1)/dx2 * (V(II+p,JJ,KK) + V(II-p,JJ,KK)) + ...
					S.w2(p+1)*S.lapc_T(2,2)/dy2 * (V(II,JJ+p,KK) + V(II,JJ-p,KK)) + ...
					S.w2(p+1)*S.lapc_T(3,3)/dz2 * (V(II,JJ,KK+p) + V(II,JJ,KK-p));
			for q = 1:S.FDn
				b(II,JJ,KK) = b(II,JJ,KK) + S.w1(p+1)*S.w1(q+1)*S.lapc_T(1,2)/dxdy * ( V(II+q,JJ+p,KK) - ...
						V(II-q,JJ+p,KK) - V(II+q,JJ-p,KK) + V(II-q,JJ-p,KK) ) + ...
						S.w1(p+1)*S.w1(q+1)*S.lapc_T(2,3)/dydz * ( V(II,JJ+q,KK+p) - ...
						V(II,JJ-q,KK+p) - V(II,JJ+q,KK-p) + V(II,JJ-q,KK-p) ) + ...
						S.w1(p+1)*S.w1(q+1)*S.lapc_T(1,3)/dzdx * ( V(II+q,JJ,KK+p) - ...
						V(II-q,JJ,KK+p) - V(II+q,JJ,KK-p) + V(II-q,JJ,KK-p) ) ;
			end
		end

	elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
		b = pseudochargeDensity_atom_cychel(V,II,JJ,KK,xin,S);
	end
end