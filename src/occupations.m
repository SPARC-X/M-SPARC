function S = occupations(S)
	
FermiEnergyEvaluator = @(lambda_f_g) fermiCalc(lambda_f_g,S);

S.lambda_f = fzero(FermiEnergyEvaluator,0);

%fprintf(2,' ------------------\n');
fprintf(' Fermi energy = %f\n',S.lambda_f);

% Calculate occupations
if S.elec_T_type == 0 % fermi-dirac smearing
	S.occ = 1./(1+exp(S.bet*(S.EigVal-S.lambda_f)));
elseif S.elec_T_type == 1 % gaussian smearing
	S.occ = 0.5 * (1.0 - erf(S.bet*(S.EigVal-S.lambda_f)));
end

end


	
function f = fermiCalc(lambda_f_g, S)
	f = 0;
	ks = 1;
	for spin = 1:S.nspin
		for kpt = 1:S.tnkpt
			if S.elec_T_type == 0 % fermi-dirac smearing
				f = f + S.occfac * sum(S.wkpt(kpt)./(1+exp(S.bet*(S.EigVal(:,ks)-lambda_f_g)))) ;
			else
				f = f + S.occfac * sum(S.wkpt(kpt)*0.5*(1.0-erf(S.bet*(S.EigVal(:,ks)-lambda_f_g)))) ;
			end
			ks = ks + 1;
		end
	end
	%f = f - S.Nelectron;
	f = f + S.NegCharge;
end