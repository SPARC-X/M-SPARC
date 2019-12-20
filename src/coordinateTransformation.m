function Y = coordinateTransformation(S, X, transfrm_typ)
% @ brief      Function to perform transformation of quantities like
%              distance, gradient etc. from one coordinate system to another
% @ authors
%         Abhiraj Sharma <asharma424@gatech.edu>
%         Qimen Xu <qimenxu@gatech.edu>
%         Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param X(Input)           N x 3 matrix
% @param Y(Output)          N x 3 matrix
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%====================================================================== 

	% Non cartesian coordinate -> cartesian coordinate
	if (strcmp(transfrm_typ,'noncart2cart_dis'))
		if S.cell_typ == 1
			Y = X;
		elseif S.cell_typ == 2
			Y = transpose(S.lat_uvec) * X';
			Y = Y';
		elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
			Y = coordinateTransformation_cychel(S,X,transfrm_typ);
		end
	elseif (strcmp(transfrm_typ,'cart2noncart_dis'))
		if S.cell_typ == 1
			Y = X;
		elseif S.cell_typ == 2
			Y = S.grad_T * X';
			Y = Y';
		elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
			Y = coordinateTransformation_cychel(S,X,transfrm_typ);
		end
	end
end