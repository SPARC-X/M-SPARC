function [S] = scf(S)
% @brief    SCF(S) implements the self consistent field calculation.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%

fprintf('\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');
fprintf(' Starting SCF iteration...\n');

outfname = S.outfname;
fileID = fopen(outfname,'a');
if (fileID == -1) 
	error('\n Cannot open file "%s"\n',outfname);
end
if(S.nspin == 1)
	fprintf(fileID,'=====================================================================\n');
	fprintf(fileID,'                    Self Consistent Field (SCF#%d)                     \n',S.Relax_iter);
	fprintf(fileID,'=====================================================================\n');
	fprintf(fileID,'Iteration     Free Energy (Ha/atom)   SCF Error        Timing (sec)\n');
else
	fprintf(fileID,'========================================================================================\n');
	fprintf(fileID,'                            Self Consistent Field (SCF#%d)                     \n',S.Relax_iter);
	fprintf(fileID,'========================================================================================\n');
	fprintf(fileID,'Iteration     Free Energy (Ha/atom)    Magnetization     SCF Error        Timing (sec)\n');
end
fclose(fileID);

% Electrostatic potential
S = poissonSolve(S, S.poisson_tol, 0);

% Exchange-correlation potential
S = exchangeCorrelationPotential(S);

% Effective potential
S.Veff = real(bsxfun(@plus,S.phi,S.Vxc));

% Initialize the mixing history vectors
S.X = zeros(S.N*S.nspin,S.MixingHistory);
S.F = zeros(S.N*S.nspin,S.MixingHistory);
S.mixing_hist_fkm1 = zeros(S.N*S.nspin,1);

if S.nspin == 1
	if S.MixingVariable == 0
		S.mixing_hist_xkm1 = S.rho;
		rho_temp = S.rho;
	else
		% for potential mixing, we store the mean-0 part only
		if S.BC == 2
			Veff_mean = mean(S.Veff);
		else
			Veff_mean = 0.0;
		end
		S.mixing_hist_xkm1 = S.Veff - Veff_mean;
		Veff_temp = S.Veff;
	end
else
	if S.MixingVariable == 0
		RHO_temp = vertcat(S.rho(:,2),S.rho(:,3));
		S.mixing_hist_xkm1 = RHO_temp;
		rho_temp = S.rho;
	else
		VEFF_temp = vertcat(S.Veff(:,1),S.Veff(:,2));
		% for potential mixing, we store the mean-0 part only
		VEFF_temp_mean = mean(VEFF_temp);
		S.mixing_hist_xkm1 = VEFF_temp - VEFF_temp_mean;
		Veff_temp = S.Veff;
	end
end

% update charges
% S.PosCharge = abs(dot(S.W, S.b));
% S.NegCharge = -S.PosCharge + S.NetCharge;

% Generate guess psi WARNING: psi is an internal matlab function
if(S.ForceCount == 1)
	rng('default'); % Initialize random number generator
	rng(1); % Specify the seed to be 1
	S.psi = rand(S.N,S.Nev,S.tnkpt*S.nspin)-0.5;
	S.upper_bound_guess_vecs = zeros(S.N,S.tnkpt*S.nspin);
	S.EigVal = zeros(S.Nev,S.tnkpt*S.nspin);
end
% S.EigVal = zeros(S.Nev,S.tnkpt*S.nspin);


% SCF LOOP 
count = 1;
count_SCF = 1;
err = 100;
max_count_first_relax = S.rhoTrigger;
max_count_gen_relax = 1;
max_scf_iter = S.MAXIT_SCF;
min_scf_iter = S.MINIT_SCF;
if max_scf_iter < min_scf_iter
	min_scf_iter = max_scf_iter;
end
S.lambda_f = 0.0;

% Spectrum bounds and filter cutoff for Chebyshev filtering
bup = zeros(S.tnkpt*S.nspin,1);
a0 = zeros(S.tnkpt*S.nspin,1);
lambda_cutoff = zeros(S.tnkpt*S.nspin,1);

% time for one scf
t_SCF = 0; 

% start scf loop
while (err > S.SCF_tol && count_SCF <= max_scf_iter || count_SCF <= min_scf_iter)
	tic_cheb = tic;
	if(count_SCF > 1)
		fprintf(' ========================= \n');
		fprintf(' Relaxation iteration: %2d \n SCF iteration number: %2d \n',S.Relax_iter,count_SCF);
		fprintf(' ========================= \n');
	else
		fprintf(' ============================================= \n');
		fprintf(' Relaxation iteration: %2d\n SCF iteration number:  1, Chebyshev cycle: %d \n',S.Relax_iter,count);
		fprintf(' ============================================= \n');
	end

	[S.upper_bound_guess_vecs,S.psi,S.EigVal,a0,bup,lambda_cutoff] = ...
	eigSolver(S,count,S.upper_bound_guess_vecs,S.psi,S.EigVal,a0,bup,lambda_cutoff);

	% Solve for Fermi energy S.lambda_f and occupations
	S = occupations(S);
	
	if (((S.ForceCount == 1) && (count >= max_count_first_relax)) ...
			|| ((S.ForceCount > 1) && (count >= max_count_gen_relax)))
		% for density mixing, can estimate energy based on input rho and
		% input veff, will recalculate energy once scf is converged
		if S.MixingVariable == 0
			[S.Etotal,S.Eband,S.Exc,S.Exc_dc,S.Eelec_dc,S.Eent] = evaluateTotalEnergy(S);
		end
		
		% Electron density
		S = electronDensity(S);

		if S.MixingVariable == 1
			% update Veff before calculating the energy %
			% Electrostatic potential
			S = poissonSolve(S, S.poisson_tol, 1);
			% Exchange-correlation potential
			S = exchangeCorrelationPotential(S);
			% Effective potential
			S.Veff = real(bsxfun(@plus,S.phi,S.Vxc));
			% Calculate energy
			[S.Etotal,S.Eband,S.Exc,S.Exc_dc,S.Eelec_dc,S.Eent] = evaluateTotalEnergy(S);
			% Calculate Self Consistency Correction to energy
			if S.nspin == 1
				S.Escc = sum((S.Veff-Veff_temp) .* S.rho .* S.W);
			else
				S.Escc = sum(sum((S.Veff-Veff_temp) .* S.rho(:,2:3),2) .* S.W);
			end
			fprintf(' Escc = %.8f\n', S.Escc);
			S.Etotal = S.Etotal + S.Escc;
        end
        fprintf(' Etot = %.8f\n', S.Etotal);
		fprintf(' Eatom = %.8f\n', S.Etotal/S.n_atm);
		
		% fprintf('\n Time for total energy calculation = %f s.',toc(t1));
		% fprintf('\n Total energy (Ha/atom): %1.9f \n',Etotal/S.n_atm) ;
		S_Debug.relax(S.Relax_iter).fermi_energy(count_SCF,1) = S.lambda_f;
		S_Debug.relax(S.Relax_iter).Etotal(count_SCF,1) = S.Etotal;
		S_Debug.relax(S.Relax_iter).Eband(count_SCF,1) = S.Eband;
		S_Debug.relax(S.Relax_iter).Exc(count_SCF,1) = S.Exc;
		S_Debug.relax(S.Relax_iter).Exc_dc(count_SCF,1) = S.Exc_dc;
		S_Debug.relax(S.Relax_iter).Eelec_dc(count_SCF,1) = S.Eelec_dc;
		S_Debug.relax(S.Relax_iter).Eent(count_SCF,1) = S.Eent;
		if S.MixingVariable == 1
			S_Debug.relax(S.Relax_iter).Escc(count_SCF,1) = S.Escc;
		end
		% Error in SCF fixed-point iteration
		if S.nspin == 1
			if S.MixingVariable == 1
				err = (norm(S.Veff - Veff_temp))/(norm(S.Veff));
			else
				err = (norm(S.rho - rho_temp))/(norm(S.rho));
			end
		else
			S.netM = dot(S.W,S.rho(:,2)) - dot(S.W,S.rho(:,3));
			fprintf('======================================\n');
			fprintf('Net magnetization in this iteration is: % .15f \n', S.netM);
			fprintf('======================================\n');
			if S.MixingVariable == 1
				VEFF = vertcat(S.Veff(:,1), S.Veff(:,2));
				err = norm(VEFF - VEFF_temp)/norm(VEFF);
			else
				RHO = vertcat(S.rho(:,2), S.rho(:,3));
				err = norm(RHO - RHO_temp)/norm(RHO);
			end
		end
		
		fprintf(' Error in SCF iteration: %.4e \n',err) ;

		% Mixing to accelerate SCF convergence
		if (err > S.SCF_tol || count_SCF < S.MINIT_SCF)
			if S.nspin == 1
				if S.MixingVariable == 1 % potential mixing
					% shift Veff and Veff_temp so they have mean 0
					if S.BC == 2
						Veff_mean = mean(S.Veff); Veff_temp_mean = mean(Veff_temp);
					else
						Veff_mean = 0; Veff_temp_mean = 0;
					end
					S.Veff = S.Veff - Veff_mean;
					Veff_temp = Veff_temp - Veff_temp_mean;
					% enter mixing
					[S,S.Veff] = mixing(S,S.Veff,Veff_temp,count_SCF);
					% shift the mean back
					S.Veff = S.Veff + Veff_mean; % the new veff for next input
					% note we add Veff_mean not Veff_temp_mean, since the one
					% used for the next input is S.Veff 
					Veff_temp = S.Veff;
				else % density mixing
					[S, S.rho] = mixing(S,S.rho,rho_temp,count_SCF);
					rho_temp = S.rho;
					% at this point rho_temp = S.rho, i.e., new input density
					% update Veff
					S = poissonSolve(S, S.poisson_tol, 1);

					% Exchange-correlation potential
					S = exchangeCorrelationPotential(S);

					% Effective potential
					S.Veff = real(bsxfun(@plus,S.phi,S.Vxc));
				end
			else
				if S.MixingVariable == 1 % potential mixing
					% shift Veff and Veff_temp so they have mean 0
					if S.BC == 2
						VEFF_mean = mean(VEFF); VEFF_temp_mean = mean(VEFF_temp);
					else
						VEFF_mean = 0; VEFF_temp_mean = 0;
					end
					VEFF = VEFF - VEFF_mean;
					VEFF_temp = VEFF_temp - VEFF_temp_mean;
					% enter mixing
					[S,VEFF] = mixing(S,VEFF,VEFF_temp,count_SCF);
					% shift the mean back
					VEFF = VEFF + VEFF_mean; % the new veff for next input
					% note we add Veff_mean not Veff_temp_mean, since the one
					% used for the next input is S.Veff
					S.Veff(:,1) = VEFF(1:S.N);
					S.Veff(:,2) = VEFF(S.N+1:end);
					Veff_temp = S.Veff;
					VEFF_temp = VEFF;
				else % density mixing
					[S, RHO] = mixing(S,RHO,RHO_temp,count_SCF);
					negrho_count = sum(RHO < 0);
					if (negrho_count > 0)
						fprintf('\nDensity got negative\n\n');
						RHO(RHO < 0) = 0;
					end
					
					S.rho(:,2) = RHO(1:S.N);
					S.rho(:,3) = RHO(S.N+1:end);
					S.rho(:,1) = S.rho(:,2) + S.rho(:,3);
					if (negrho_count > 0)
						scal = (-S.NegCharge/dot(S.W,S.rho(:,1)));
						S.rho = scal*S.rho;
					end
					rho_temp = S.rho;
					RHO_temp = vertcat(rho_temp(:,2),rho_temp(:,3));
					% at this point rho_temp = S.rho, i.e., new input density
					% update Veff
					S = poissonSolve(S, S.poisson_tol, 1);

					% Exchange-correlation potential
					S = exchangeCorrelationPotential(S);

					% Effective potential
					S.Veff = real(bsxfun(@plus,S.phi,S.Vxc));
				end
			end
		end
		
		% write to output file
		fileID = fopen(S.outfname,'a');
		if (fileID == -1) 
			error('\n Cannot open file "%s"\n',S.outfname);
		end
		if t_SCF == 0
			t_SCF = toc(tic_cheb);
		end
		
		if S.nspin == 1
			fprintf(fileID,'%-6d      %18.10E        %.3E        %.3f\n', ...
					count_SCF, S.Etotal/S.n_atm, err, t_SCF);
		else
			fprintf(fileID,'%-6d      %18.10E        %11.4E        %.3E        %.3f\n', ...
					count_SCF, S.Etotal/S.n_atm, S.netM, err, t_SCF);
		end
		t_SCF = 0; % reset SCF timer
		fclose(fileID);
		count_SCF = count_SCF + 1 ;
	end
	scf_runtime = toc(tic_cheb);
	t_SCF = t_SCF + scf_runtime; % add chebyshev filtering time to SCF time
	S_Debug.relax(S.Relax_iter).scf_runtime(count_SCF,1) = scf_runtime;
	fprintf(' This SCF iteration took %.3f s.\n\n', scf_runtime);
	count = count + 1;
end

S_Debug.relax(S.Relax_iter).scf_flag = 0; % scf_flag being 0 means it's converged
if (count_SCF == max_scf_iter + 1)
	disp(' SCF did not converge. Maximum iterations reached!')
	S_Debug.relax(S.Relax_iter).scf_flag = 1; % scf_flag being 1 means it's not converged
end

S_Debug.relax(S.Relax_iter).count_SCF = count_SCF - 1;
fprintf('\n Finished SCF iteration in %d steps!\n', (count_SCF - 1));
fprintf(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');

% write to output file
fileID = fopen(S.outfname,'a');
if (fileID == -1) 
	error('\n Cannot open file "%s"\n',S.outfname);
end
fprintf(fileID,'Total number of SCF: %-6d\n',count_SCF-1);
fclose(fileID);

end

