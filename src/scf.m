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
if S.spin_typ == 0
	fprintf(fileID,'=====================================================================\n');
	fprintf(fileID,'                    Self Consistent Field (SCF#%d)                     \n',S.Relax_iter);
	fprintf(fileID,'=====================================================================\n');
	fprintf(fileID,'Iteration     Free Energy (Ha/atom)   SCF Error        Timing (sec)\n');
elseif S.spin_typ == 1
	fprintf(fileID,'========================================================================================\n');
	fprintf(fileID,'                            Self Consistent Field (SCF#%d)                     \n',S.Relax_iter);
	fprintf(fileID,'========================================================================================\n');
	fprintf(fileID,'Iteration     Free Energy (Ha/atom)    Magnetization     SCF Error        Timing (sec)\n');
elseif S.spin_typ == 2
    fprintf(fileID,'======================================================================================================================\n');
	fprintf(fileID,'                                              Self Consistent Field (SCF#%d)                                          \n',S.Relax_iter);
	fprintf(fileID,'======================================================================================================================\n');
	fprintf(fileID,'Iteration     Free Energy (Ha/atom)            Magnetization (x,y,z,tot)                 SCF Error        Timing (sec)\n');
end
fclose(fileID);

% Electrostatic potential
S = poissonSolve(S, S.poisson_tol, 0);

% Exchange-correlation potential
S.countPotential = -1;
S = exchangeCorrelationPotential(S);

% Effective potential
S.Veff = real(bsxfun(@plus,S.phi,S.Vxc));
if S.spin_typ == 2
    S.Veff_nc = S.Vxc_nc;
    S.Veff_nc(:,1:2) = S.Veff_nc(:,1:2) + S.phi;
end

% Initialize the mixing history vectors
% TODO: change it for non-collinear
S.X = zeros(S.N*S.nspden,S.MixingHistory);
S.F = zeros(S.N*S.nspden,S.MixingHistory);
S.mixing_hist_fkm1 = zeros(S.N*S.nspden,1);

if S.spin_typ == 0
    if S.MixingVariable == 0
		S.mixing_hist_xkm1 = S.rho;		
	else
		% for potential mixing, we store the mean-0 part only
		[Veff_shift] = shifting_Veff(S.Veff,S.BC,S.spin_typ);
		S.mixing_hist_xkm1 = Veff_shift;
    end
elseif S.spin_typ == 1
    if S.MixingVariable == 0		
		S.mixing_hist_xkm1 = reshape(S.rho(:,2:3),[],1);
    else
        % for potential mixing, we store the mean-0 part only
        [Veff_shift] = shifting_Veff(S.Veff,S.BC,S.spin_typ);
		S.mixing_hist_xkm1 = Veff_shift(:);
    end
elseif S.spin_typ == 2
    if S.MixingVariable == 0
        % [rho_tot,mx,my,mz]
		S.mixing_hist_xkm1 = [S.rho(:,1); reshape(S.mag(:,1:3),[],1)];
    else        
        % for potential mixing, we store the mean-0 part only
        [Veff_shift] = shifting_Veff(S.Veff_nc,S.BC,S.spin_typ);
		S.mixing_hist_xkm1 = Veff_shift(:);
    end
end

% Generate guess psi WARNING: psi is an internal matlab function
if S.ForceCount == 1
	rng('default'); % Initialize random number generator
	rng(1); % Specify the seed to be 1
	S.psi = rand(S.N*S.nspinor,S.Nev,S.tnkpt)-0.5;
	S.upper_bound_guess_vecs = zeros(S.N*S.nspinor_eig,S.num_eig);
	S.EigVal = zeros(S.Nev*S.nspin,S.tnkpt);
end

if S.usefock == 1
    scf_tol_init = S.SCF_tol_init;
else
    scf_tol_init = S.SCF_tol;
end


S.rhoin = [];
S.rhoout = [];
S.rhomix = [];

S.lambda_f = 0.0;
S = scf_loop(S,scf_tol_init);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%          SCF LOOP Exact Exchange            %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exact exchange potential 
if S.usefock == 1
    S.usefock = S.usefock+1;
else
    return;
end

% Exchange-correlation potential
S = exchangeCorrelationPotential(S);

% Effective potential
S.Veff = real(bsxfun(@plus,S.phi,S.Vxc));
if S.spin_typ == 2
    S.Veff_nc = S.Vxc_nc;
    S.Veff_nc(:,1:2) = S.Veff_nc(:,1:2) + S.phi;
end

% Note: No need to update mixing_hist_xk. Done in mixing of first iteration.

% Exact exchange potential parameters
err_Exx = S.FOCK_TOL+1;
count_Exx = 1;

S.lambda_f = 0.0;
while(count_Exx <= S.MAXIT_FOCK)
    % Store orbitals and occupations for outer loop
    S.psi_outer = S.psi;                                                                                                                                                                                                                                                                                                 
    S.occ_outer = S.occ;
    
    if S.ACEFlag == 1
        S = ace_operator(S);
    end
    
    % Calculate estimation of Exact Exchange energy
    S = evaluateExactExchangeEnergy(S);
    
    Eexx_pre = S.Eex;
    
    fileID = fopen(S.outfname,'a');
    fprintf(fileID, '\nNo.%d Exx outer loop:\n', count_Exx);
    fclose(fileID);
    
    % Start SCF with hybrid functional
    S = scf_loop(S,S.SCF_tol,count_Exx);
    
    % update exact exchange energy
    S.Etotal = S.Etotal + 2*S.Eex;
    S.Exc = S.Exc - S.Eex;
    % Calculate accurate Exact Exchange energy
    S = evaluateExactExchangeEnergy(S);
    % update exact exchange energy
    S.Etotal = S.Etotal - 2*S.Eex;
    S.Exc = S.Exc + S.Eex;
    
    err_Exx = abs(S.Eex - Eexx_pre)/S.n_atm;        
    fprintf(' Exx outer loop error: %.4e \n',err_Exx) ;
    fileID = fopen(S.outfname,'a');
    fprintf(fileID, 'Exx outer loop error: %.4e \n', err_Exx);
    fclose(fileID);

    if err_Exx < S.FOCK_TOL && count_Exx >= S.MINIT_FOCK
        break;
    end

    count_Exx = count_Exx + 1;
end % end of Vxx loop

fprintf('\n Finished outer loop in %d steps!\n', (count_Exx - 1));
fprintf(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');

if count_Exx > S.MAXIT_FOCK && err_Exx > S.FOCK_TOL
    disp(' Exact Exchange outer loop did not converge. Maximum iterations reached!');
    fileID = fopen(S.outfname,'a');
    fprintf(fileID, ...
        ' Warning: Exact Exchange outer loop did not converge. Maximum iterations reached!\n ');
    fclose(fileID);
end

% make sure next scf starts with normal scf
S.usefock = S.usefock+1;
end


function S = scf_loop(varargin)
S = varargin{1};
scf_tol = varargin{2};
count_Exx = -1;

if nargin == 3
	count_Exx = varargin{3};
elseif nargin > 3 || nargin < 2
	error('Too many input arguments.');
end

% SCF LOOP 
count = 1;
count_SCF = 1;
if S.ForceCount > 1 || S.usefock > 1    
    S.rhoTrigger = 1;
end
if count_Exx > 0
    S.MINIT_SCF = 1;
end

% Spectrum bounds and filter cutoff for Chebyshev filtering
bup = zeros(S.num_eig,1);
a0 = zeros(S.num_eig,1);
lambda_cutoff = zeros(S.num_eig,1);

% time for one scf
t_SCF = 0; 
if S.spin_typ == 0
	if S.MixingVariable == 0
		rho_temp = S.rho;
    else
		Veff_temp = S.Veff;
	end
elseif S.spin_typ == 1
    if S.MixingVariable == 0
        rho_temp = S.rho(:,2:3);
    else
        Veff_temp = S.Veff;
    end
elseif S.spin_typ == 2
    if S.MixingVariable == 0
        rho_temp = S.rho(:,2:3); 
        rho_nc_temp = [S.rho(:,1) S.mag(:,1:3)]; 
    else
        Veff_temp = S.Veff;
        Veff_nc_temp = S.Veff_nc;
    end
end

% start scf loop
while count_SCF <= S.MAXIT_SCF
	tic_cheb = tic;
    if count_SCF == 1
        Ncheb = S.rhoTrigger;
    else
        Ncheb = 1;
    end
    
    for ncheb = 1:Ncheb
        if count_SCF > 1
            fprintf(' ========================= \n');
            if S.usefock > 1
                fprintf(' Outer loop iteration number: %2d\n', count_Exx);
            end
            fprintf(' Relaxation iteration: %2d \n SCF iteration number: %2d \n',S.Relax_iter,count_SCF);
            fprintf(' ========================= \n');
        else
            fprintf(' ============================================= \n');
            if S.usefock > 1
                fprintf(' Outer loop iteration number: %2d\n', count_Exx);
            end
            fprintf(' Relaxation iteration: %2d\n SCF iteration number:  1, Chebyshev cycle: %d \n',S.Relax_iter,count);
            fprintf(' ============================================= \n');
        end
    
        [S.upper_bound_guess_vecs,S.psi,S.EigVal,a0,bup,lambda_cutoff] = ...
        eigSolver(S,count,S.upper_bound_guess_vecs,S.psi,S.EigVal,a0,bup,lambda_cutoff);

        % Solve for Fermi energy S.lambda_f and occupations
        S = occupations(S);
        count = count + 1;
    end
    
    % for density mixing, us rho_in to estimate total energy
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
        if S.spin_typ == 2
            S.Veff_nc = S.Vxc_nc;
            S.Veff_nc(:,1:2) = S.Veff_nc(:,1:2) + S.phi;
        end
        
        % Calculate energy
        [S.Etotal,S.Eband,S.Exc,S.Exc_dc,S.Eelec_dc,S.Eent] = evaluateTotalEnergy(S);
        % Calculate Self Consistency Correction to energy
        if S.spin_typ == 0
            S.Escc = sum((S.Veff-Veff_temp) .* S.rho .* S.W);
        else
            S.Escc = sum(sum((S.Veff-Veff_temp).* S.rho(:,2:3),2) .* S.W);
        end
        fprintf(' Escc = %.8f\n', S.Escc);
        S.Etotal = S.Etotal + S.Escc;
    end
    fprintf(' Etot = %.8f\n', S.Etotal);
    fprintf(' Eatom = %.8f\n', S.Etotal/S.n_atm);

    %%%%%%%%%%%%%%%%%%%%%%%%%   debug info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %%%%%%%%%%%%%%%%%%%%%%%%%   debug info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Error in SCF fixed-point iteration
    if S.spin_typ == 0
        if S.MixingVariable == 1
            err = (norm(S.Veff - Veff_temp))/(norm(S.Veff));
        else
            err = (norm(S.rho - rho_temp))/(norm(S.rho));
        end
    else
        S.netM = sum(S.mag)*S.dV;
        fprintf('======================================\n');
        fprintf('Net magnetization in this iteration is: ');
        fprintf('%.6f ', S.netM);
        fprintf('\n');
        fprintf('======================================\n');
        if S.MixingVariable == 1
            Veff = S.Veff;
            err = norm(reshape(Veff - Veff_temp,[],1))/norm(Veff(:));
        else
            rho = S.rho(:,2:3);
            err = norm(reshape(rho - rho_temp,[],1))/norm(rho(:));
        end
    end

    fprintf(' Error in SCF iteration: %.4e \n',err) ;

    % Mixing to accelerate SCF convergence    
    if S.spin_typ == 0
        if S.MixingVariable == 1 % potential mixing
            % shift Veff and Veff_temp so they have mean 0            
            [S.Veff,Veff_mean] = shifting_Veff(S.Veff,S.BC,S.spin_typ);
            [Veff_temp] = shifting_Veff(Veff_temp,S.BC,S.spin_typ);
            
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
        end
    elseif S.spin_typ == 1
        if S.MixingVariable == 1 % potential mixing
            [Veff,Veff_mean] = shifting_Veff(Veff,S.BC,S.spin_typ);
            [Veff_temp] = shifting_Veff(Veff_temp,S.BC,S.spin_typ);

            % enter mixing
            [S,Veff(:)] = mixing(S,Veff(:),Veff_temp(:),count_SCF);

            Veff = Veff + Veff_mean; % the new veff for next input
            S.Veff = Veff;
            Veff_temp = Veff;
        else
            [S, rho(:)] = mixing(S,rho(:),rho_temp(:),count_SCF);
            
            S.rho(:,2:3) = rho;
            S.rho(:,1) = S.rho(:,2) + S.rho(:,3);
            S.mag = S.rho(:,2) - S.rho(:,3);
            rho_temp = S.rho(:,2:3);
        end
    else
        if S.MixingVariable == 1 % potential mixing
            Veff_nc = S.Veff_nc;
            [Veff_nc,Veff_mean] = shifting_Veff(Veff_nc,S.BC,S.spin_typ);
            [Veff_nc_temp] = shifting_Veff(Veff_nc_temp,S.BC,S.spin_typ);

            % enter mixing
            [S,Veff_nc(:)] = mixing(S,Veff_nc(:),Veff_nc_temp(:),count_SCF);

            % shift the mean back
            Veff_nc = Veff_nc + Veff_mean; % the new veff for next input
            S.Veff_nc = Veff_nc;

            Veff_nc_temp = Veff_nc;
            Veff_temp = S.Veff;
        else
            rho_nc = [S.rho(:,1) S.mag(:,1:3)];
            [S, rho_nc(:)] = mixing(S,rho_nc(:),rho_nc_temp(:),count_SCF);
            rho_nc_temp = rho_nc;
            
            S.mag(:,1:3) = rho_nc(:,2:4);
            S.mag(:,4) = sqrt(sum(S.mag(:,1:3).*S.mag(:,1:3),2));
            S.rho = [rho_nc(:,1) 0.5*(rho_nc(:,1)+S.mag(:,4)) 0.5*(rho_nc(:,1)-S.mag(:,4))];
            rho_temp = S.rho(:,2:3);
        end
    end    
    
    if S.MixingVariable == 0
        % update Veff
        S = poissonSolve(S, S.poisson_tol, 1);

        % Exchange-correlation potential
        S = exchangeCorrelationPotential(S);

        % Effective potential            
        S.Veff = real(bsxfun(@plus,S.phi,S.Vxc));
        if S.spin_typ == 2
            S.Veff_nc = S.Vxc_nc;
            S.Veff_nc(:,1:2) = S.Veff_nc(:,1:2) + S.phi;
        end
    end

    % write to output file
    fileID = fopen(S.outfname,'a');
    if (fileID == -1) 
        error('\n Cannot open file "%s"\n',S.outfname);
    end
    
    scf_runtime = toc(tic_cheb);
	t_SCF = t_SCF + scf_runtime; % add chebyshev filtering time to SCF time
	S_Debug.relax(S.Relax_iter).scf_runtime(count_SCF,1) = scf_runtime;
	fprintf(' This SCF iteration took %.3f s.\n\n', scf_runtime);	

    if S.spin_typ == 0
        fprintf(fileID,'%-6d      %18.10E        %.3E        %.3f\n', ...
                count_SCF, S.Etotal/S.n_atm, err, scf_runtime);
    elseif S.spin_typ == 1
        fprintf(fileID,'%-6d      %18.10E        %11.4E        %.3E        %.3f\n', ...
                count_SCF, S.Etotal/S.n_atm, S.netM, err, scf_runtime);
    elseif S.spin_typ == 2
        fprintf(fileID,'%-6d      %18.10E    %11.4E, %11.4E, %11.4E, %11.4E     %.3E        %.3f\n', ...
                count_SCF, S.Etotal/S.n_atm, S.netM, err, scf_runtime);
    end
    fclose(fileID);
    count_SCF = count_SCF + 1 ;
    
    if err < scf_tol && count_SCF > S.MINIT_SCF
        break;
    end
end

S_Debug.relax(S.Relax_iter).scf_flag = 0; % scf_flag being 0 means it's converged
if (count_SCF == S.MAXIT_SCF + 1)
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



function [Veff_shift,Veff_mean] = shifting_Veff(Veff,BC,spin_typ)
ncol = size(Veff,2);
if spin_typ ~= 2
    if BC == 2
        Veff_mean = mean(Veff(:));
    else
        Veff_mean = 0;
    end
    Veff_shift = reshape(Veff(:) - Veff_mean,[],ncol);
else
    if BC == 2
        Veff_mean = mean(Veff);
    else
        Veff_mean = zeros(1,4);
    end
    Veff_shift = Veff - Veff_mean;
end        
end

