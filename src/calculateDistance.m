function dd = calculateDistance(X,Y,Z,X_ref,Y_ref,Z_ref,S)
% @ brief      Function to perform distance calculation
% @ authors
%         Abhiraj Sharma <asharma424@gatech.edu>
%         Qimen Xu <qimenxu@gatech.edu>
%         Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param X          X coordinate of the nodes
% @param Y          Y coordinate of the nodes
% @param Z          Z coordinate of the nodes
% @param X_ref       x-coordinate of the reference point
% @param Y_ref       y-coordinate of the reference point
% @param Z_ref       z-coordinate of the reference point
% @param dd          distance of the nodes from the reference point
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%====================================================================== 	
	if S.cell_typ == 1
		XX = bsxfun(@minus,X,X_ref);
		YY = bsxfun(@minus,Y,Y_ref);
		ZZ = bsxfun(@minus,Z,Z_ref);
		dd = sqrt(XX.^2 + YY.^2 + ZZ.^2);
	elseif S.cell_typ == 2
		XX = bsxfun(@minus,X,X_ref);
		YY = bsxfun(@minus,Y,Y_ref);
		ZZ = bsxfun(@minus,Z,Z_ref);
		dd = sqrt(S.metric_T(1,1)*XX.^2 + S.metric_T(1,2)*(XX.*YY) + S.metric_T(1,3)*(XX.*ZZ) + ...
				  S.metric_T(2,2)*YY.^2 + S.metric_T(2,3)*(YY.*ZZ) + S.metric_T(3,3)*ZZ.^2);
	elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
		dd = calculateDistance_cychel(X,Y,Z,X_ref,Y_ref,Z_ref,S);
	end
end