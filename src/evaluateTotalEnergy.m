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

if S.nspin == 1
	if S.xc == 0 % LDA_PW
		C2 = 0.73855876638202 ; % constant for exchange energy
		%rho = rho+(1e-50) ; % to avoid divide by zero error
		p = 1 ;
		A = 0.031091 ;
		alpha1 = 0.21370 ;
		beta1 = 7.5957 ;
		beta2 = 3.5876 ;
		beta3 = 1.6382 ;
		beta4 = 0.49294 ;
		CEnergyPotential = (0.75./(pi*rho)).^(1/3) ;
		CEnergyPotential = -2*A*(1+alpha1*CEnergyPotential).*log(1+1./(2*A*( beta1*(CEnergyPotential.^0.5) ...
		   + beta2*CEnergyPotential + beta3*(CEnergyPotential.^1.5) + beta4*(CEnergyPotential.^(p+1.0))))) ;
		%rho = rho-(1e-50) ;
		Exc = sum(CEnergyPotential.*rho.*S.W) - C2*sum((rho.^(4/3)).*S.W) ;
	elseif S.xc == 1 % LDA_PZ
		A = 0.0311;
		B = -0.048 ;
		C = 0.002 ;
		D = -0.0116 ;
		gamma1 = -0.1423 ;
		beta1 = 1.0529 ;
		beta2 = 0.3334 ;
		C2 = 0.73855876638202;
		%rho = rho+(1e-50) ; % to avoid divide by zero error
		CEnergyPotential = (0.75./(pi*rho)).^(1/3) ;
		islt1 = (CEnergyPotential < 1.0);
		CEnergyPotential(islt1) = A * log(CEnergyPotential(islt1)) + B ...
		   + C * CEnergyPotential(islt1) .* log(CEnergyPotential(islt1)) ...
		   + D * CEnergyPotential(islt1);
		CEnergyPotential(~islt1) = gamma1./(1.0+beta1*sqrt(CEnergyPotential(~islt1))+beta2*CEnergyPotential(~islt1));
		%rho = rho-(1e-50) ;
		Exc = sum(CEnergyPotential.*rho.*S.W) - C2*sum((rho.^(4/3)).*S.W) ;
	elseif S.xc == 2
		Exc = sum(S.e_xc.*rho.*S.W);
	end
	% Exchange-correlation energy double counting correction
	Exc_dc = sum(S.Vxc.*rho.*S.W) ;
else
	Exc = sum(S.e_xc.*rho(:,1).*S.W);
	% Exchange-correlation energy double counting correction
	Exc_dc = sum(sum(S.Vxc.*rho(:,2:3),2).*S.W) ;
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
Etot = Eband + Exc - Exc_dc + Eelec_dc - S.Eself + S.E_corr + Eent;

%fprintf(2,' ------------------\n');
% fprintf(' Eband = %.8f\n', Eband);
% fprintf(' Exc = %.8f\n', Exc);
% fprintf(' Exc_dc = %.8f\n', Exc_dc);
% fprintf(' Eelec_dc = %.8f\n', Eelec_dc);
% fprintf(' Eent = %.8f\n', Eent);
% fprintf(' E_corr = %.8f\n', S.E_corr);
% fprintf(' Eself = %.8f\n', S.Eself);
% fprintf(' Etot = %.8f\n', Etot);
%fprintf(2,' ------------------\n');

end
