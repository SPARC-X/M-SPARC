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
	fprintf(fileID,'Iteration     Free Energy (Ha/atom)            Magnetization (tot,x,y,z)                 SCF Error        Timing (sec)\n');
end
fclose(fileID);

% Electrostatic potential
S = poissonSolve(S, S.poisson_tol, 0);

% Exchange-correlation potential
S.countPotential = -1;
S = exchangeCorrelationPotential(S);

% Effective potential
S = calculate_effective_potential(S);

% Initialize the mixing history vectors
S.X = zeros(S.N*S.nspden,S.MixingHistory);
S.F = zeros(S.N*S.nspden,S.MixingHistory);
S.mixing_hist_xkm1 = zeros(S.N*S.nspden,1);
S.mixing_hist_fkm1 = zeros(S.N*S.nspden,1);

% initialize history
if S.MixingVariable == 0
    S.mixing_hist_xkm1(1:S.N) = S.rho(:,1);
    if S.spin_typ == 1
        S.mixing_hist_xkm1(S.N+1:end) = S.mag;
    end
    if S.spin_typ == 2
        S.mixing_hist_xkm1(S.N+1:end) = S.mag(1+S.N:end)';
    end
    S.mixing_hist_xk = S.mixing_hist_xkm1;
else 
    S.mixing_hist_xkm1(:) = S.Veff(:);    
    S.mixing_hist_xkm1 = shifting_Veff(S,S.mixing_hist_xkm1,S.nspden,[],0,-1);                        
    S.mixing_hist_xk = S.mixing_hist_xkm1;
end


% Generate guess psi WARNING: psi is an internal matlab function
if S.ForceCount == 1
	rng('default'); % Initialize random number generator
	rng(1); % Specify the seed to be 1
	S.psi = rand(S.N*S.nspinor,S.Nev,S.tnkpt)-0.5;
	S.upper_bound_guess_vecs = zeros(S.N*S.nspinor_eig,S.tnkpt*S.nspin);
	S.EigVal = zeros(S.Nev,S.tnkpt*S.nspin);
end

if S.usefock == 1
    scf_tol_init = S.SCF_tol_init;
else
    scf_tol_init = S.SCF_tol;
end

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
S = calculate_effective_potential(S);

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
bup = zeros(S.tnkpt*S.nspin,1);
a0 = zeros(S.tnkpt*S.nspin,1);
lambda_cutoff = zeros(S.tnkpt*S.nspin,1);

% time for one scf
t_SCF = 0; 
% start scf loop
while count_SCF <= S.MAXIT_SCF
	tic_cheb = tic;
    if count_SCF == 1 && S.ForceCount == 1 && S.usefock <= 1
        Nchefsi = S.rhoTrigger;
    else
        Nchefsi = S.nchefsi;
    end
    
    if S.MixingVariable == 0
        if S.spin_typ == 0
            rho_in = S.rho;
        else
            rho_in = S.rho(:,2:3);
        end
    else 
        if S.spin_typ == 2
            Veff_in = S.Veff_dia;
        else
            Veff_in = S.Veff;
        end
    end
    
    for nchefsi = 1:Nchefsi
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
        S = calculate_effective_potential(S);
        
        % Calculate energy
        [S.Etotal,S.Eband,S.Exc,S.Exc_dc,S.Eelec_dc,S.Eent] = evaluateTotalEnergy(S);
        % Calculate Self Consistency Correction to energy
        if S.spin_typ == 0
            S.Escc = sum((S.Veff-Veff_in) .* S.rho .* S.W);
        elseif S.spin_typ == 1
            S.Escc = sum(sum((S.Veff-Veff_in).* S.rho(:,2:3),2) .* S.W);
        else
            S.Escc = sum(sum((S.Veff_dia-Veff_in).* S.rho(:,2:3),2) .* S.W);
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
    if S.MixingVariable == 0
        if S.spin_typ == 0
            err = (norm(S.rho - rho_in))/(norm(S.rho));
        else
            err = norm(reshape(S.rho(:,2:3) - rho_in,[],1))/norm(reshape(S.rho(:,2:3),[],1));
        end
    else
        if S.spin_typ == 2
            err = norm(reshape(S.Veff_dia - Veff_in,[],1))/norm(S.Veff_dia(:));
        else
            err = norm(reshape(S.Veff - Veff_in,[],1))/norm(S.Veff(:));
        end
    end
    
    if S.spin_typ > 0
        S.netM = sum(S.mag)*S.dV;
        fprintf('======================================\n');
        fprintf('Net magnetization in this iteration is: ');
        fprintf('%.6f ', S.netM);
        fprintf('\n');
        fprintf('======================================\n');
    end

    fprintf(' Error in SCF iteration: %.4e \n',err) ;

    % Mixing to accelerate SCF convergence    
    if S.MixingVariable == 0
        % get rho_out
        rho_out = zeros(S.N*S.nspden,1);
        rho_out(1:S.N) = S.rho(:,1);
        if S.spin_typ == 1
            rho_out(1+S.N:end) = S.mag;
        end
        if S.spin_typ == 2
            rho_out(1+S.N:end) = S.mag(S.N+1:end)';
        end
        % mixing
        [S, rho_out] = mixing(S,rho_out,S.mixing_hist_xk,count_SCF);
        % update
        S.rho(:,1) = rho_out(1:S.N);
        if S.spin_typ == 1
            S.mag = rho_out(S.N+1:end);
            S.rho(:,2) = 0.5*(S.rho(:,1)+S.mag);
            S.rho(:,3) = 0.5*(S.rho(:,1)-S.mag);
        end
        if S.spin_typ == 2
            S.mag(:,2:4) = reshape(rho_out(S.N+1:end),[],3);
            S.mag(:,1) = sqrt(S.mag(:,2).^2+S.mag(:,3).^2+S.mag(:,4).^2);
            S.rho(:,2) = 0.5*(S.rho(:,1)+S.mag(:,1));
            S.rho(:,3) = 0.5*(S.rho(:,1)-S.mag(:,1));
        end
    else
        Veff_out = S.Veff(:);        
        % move mean
        [Veff_out,Veff_mean] = shifting_Veff(S,Veff_out,S.nspden,[],0,-1);     
        % mixing
        [S,Veff_out] = mixing(S,Veff_out,S.mixing_hist_xk,count_SCF);
        % shift back
        Veff_out = shifting_Veff(S,Veff_out,S.nspden,Veff_mean,1,1);  
        % update
        S.Veff = reshape(Veff_out,[],S.nspden);
    end
    
    if S.MixingVariable == 0
        % update Veff
        S = poissonSolve(S, S.poisson_tol, 1);
        % Exchange-correlation potential
        S = exchangeCorrelationPotential(S);
        % Effective potential            
        S = calculate_effective_potential(S);
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



function [Veff,Veff_mean] = shifting_Veff(S,Veff,ncol,Veff_mean,option,dir)
if S.BC ~= 2
    Veff_mean = zeros(ncol,1);
    return;
end

Veff = reshape(Veff,[],ncol);
if option == 0
    Veff_mean = mean(Veff);
    if S.spin_typ == 1
        Veff_mean = mean(Veff_mean);
    end
end

Veff = Veff + dir*Veff_mean;
Veff = Veff(:);
end



function [S] = calculate_effective_potential(S)
if S.spin_typ == 2
    S.Veff = S.Vxc_nc;
    S.Veff(:,1:2) = S.Veff(:,1:2) + S.phi;
    S.Veff_dia = real(S.Vxc + S.phi);
else
    S.Veff = S.Vxc + S.phi;
end
S.Veff = real(S.Veff);
end
