function [Etot,Eband,Exc,Exc_dc,Eelec_dc,Eent] = evaluateTotalEnergy(S)
% @brief    EVALUATETOTALENERGY calculates the total energy based on the
%           output density at each SCF.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%

% Band structure energy
Eband = 0;
ks = 1;
for spin = 1:S.nspin
	for kpt = 1:S.tnkpt
		Eband = Eband + S.occfac * S.wkpt(kpt) * sum(S.EigVal(:,ks).*S.occ(:,ks)) ;
		ks = ks + 1;
	end
end

% Exchange-correlation energy
rho = S.rho;
% Check if density is too small
INDX_zerorho = (rho < S.xc_rhotol);
rho(INDX_zerorho) = S.xc_rhotol;

if S.spin_typ == 0
    Exc = sum(S.e_xc.*(rho+S.rho_Tilde_at).*S.W);
    if (S.vdWDFFlag == 1) || (S.vdWDFFlag == 2) % add vdW energy in Exc
        Exc = Exc + S.vdWenergy; 
    end
	% Exchange-correlation energy double counting correction
	Exc_dc = sum(S.Vxc.*rho.*S.W) ;
    if (S.countPotential > 0) && (S.xc == 4) % S.xc == 4 SCAN functional
        Eext_scan_dc = sum(S.VxcScan3.*S.tau.*S.W);
        Exc_dc = Exc_dc + Eext_scan_dc;
    end
else
	Exc = sum(S.e_xc.*(rho(:,1)+S.rho_Tilde_at).*S.W);
	% Exchange-correlation energy double counting correction
	Exc_dc = sum(sum(S.Vxc.*rho(:,2:3),2).*S.W);
    if (S.vdWDFFlag == 1) || (S.vdWDFFlag == 2) % add vdW energy in Exc
		Exc = Exc + S.vdWenergy; 
    end
    if (S.countPotential > 0) && (S.xc == 4) % S.xc == 4 SCAN functional
        Eext_scan_dc = sum(sum(S.VxcScan3.*S.tau(:, 2:3).*S.W));
        Exc_dc = Exc_dc + Eext_scan_dc;
    end
end

% Electrostatic energy double counting correction
Eelec_dc = 0.5*sum((S.b-S.rho(:,1)).*S.phi.*S.W);

% Electronic entropy
Eent = 0 ;
ks = 1;
for spin = 1:S.nspin
	for kpt = 1:S.tnkpt
		if S.elec_T_type == 0 % fermi-dirac smearing
			Eent_v = S.occfac*(1/S.bet)*(S.occ(:,ks).*log(S.occ(:,ks))+(1-S.occ(:,ks)).*log(1-S.occ(:,ks)));
			Eent_v(isnan(Eent_v)) = 0.0 ;
		elseif S.elec_T_type == 1 % gaussian smearing
			Eent_v = -S.occfac*(1/S.bet)*1/(2*sqrt(pi)) .* exp(-(S.bet * (S.EigVal(:,ks)-S.lambda_f)).^2);
		end
		Eent = Eent + S.wkpt(kpt)*sum(Eent_v);
		ks = ks + 1;
	end
end

% Total free energy
if S.usefock < 2
    % Total free energy
    Etot = Eband + Exc - Exc_dc + Eelec_dc - S.Eself + S.E_corr + Eent;
else
    Exc = Exc + S.Eex;
    % Total free energy
    Etot = Eband + Exc - Exc_dc + Eelec_dc - S.Eself + S.E_corr + Eent - 2*S.Eex;
end


fprintf(2,' ------------------\n');
fprintf(' Eband = %.8f\n', Eband);
fprintf(' Exc = %.8f\n', Exc);
fprintf(' Exc_dc = %.8f\n', Exc_dc);
fprintf(' Eelec_dc = %.8f\n', Eelec_dc);
fprintf(' Eent = %.8f\n', Eent);
fprintf(' E_corr = %.8f\n', S.E_corr);
fprintf(' Eself = %.8f\n', S.Eself);
if S.usefock > 1
    fprintf(' Eex = %.8f\n', S.Eex);
end
fprintf(' Etot = %.8f\n', Etot);
fprintf(2,' ------------------\n');

end
