function S = electronDensity(S)
% @brief    Calculate new density based on the new states.
%
% @param S  A struct that contains the relevant fields.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%
S.rho = 0 * S.rho;
if S.nspin == 1
	for kpt =1:S.tnkpt
		S.rho = S.rho + 2*S.wkpt(kpt)* sum( bsxfun(@times,S.psi(:,:,kpt).*conj(S.psi(:,:,kpt)),S.occ(:,kpt)'),2);
	end
	S.rho = real(S.rho);
else
	ks = 1;
	for spin =1:S.nspin
		for kpt =1:S.tnkpt
			S.rho(:,spin+1) = S.rho(:,spin+1) + S.wkpt(kpt)* sum( bsxfun(@times,S.psi(:,:,ks).*conj(S.psi(:,:,ks)),S.occ(:,ks)'),2);
			ks = ks + 1;
		end
	end
	S.rho(:,2:3) = real(S.rho(:,2:3));
	S.rho(:,1) = S.rho(:,2) + S.rho(:,3);
end
	
end