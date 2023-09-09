function [upper_bound_guess_vecs,psi,EigVal,a0,bup,lambda_cutoff] = eigSolver(S,count,upper_bound_guess_vecs,psi,EigVal,a0,bup,lambda_cutoff)
% @brief    Solve the linearized Hamiltonian eigenproblem.
%
% @param count                  : SCF iteration count.
% @param upper_bound_guess_vecs : guess vector for maximum eigenvector of
%                                 the linearized Hamiltonian.
% @param psi                    : eigenvectors of the previous linearized
%                                 Hamiltonian.
% @param EigVal                 : previous eigenvalues
% @param a0                     : lower bound for Chebyshev filtering (1st SCF).
% @param bup                    : upper bound for Chebyshev filtering (1st SCF).
% @param lambda_cutoff          : cut-off for chebyshev filtering (1st SCF).
%
% @authors	Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%

% eigen problem size
N = S.N*S.nspinor_eig; 
% reorganize psi, EigVal to make each chunk of data for one eigen problem
% while S.psi remains original format
if S.spin_typ == 1
    psi = cat(3,psi(1:S.N,:,:),psi(1+S.N:2*S.N,:,:));
end

if S.parallel ~= 1
    for ks = 1:S.tnkpt*S.nspin
        if ks <= S.tnkpt
            kpt = ks;
            spin = 1;
        else
            kpt = ks - S.tnkpt;
            spin = 2;
        end
        % get Veff 
        if S.spin_typ == 2
            Veff = S.Veff;
        else
            Veff = S.Veff(:,spin);
        end

        rng('default'); % Initialize random number generator
		rng(ks+1);
		opts = struct('maxit', 100, 'tol', S.TOL_LANCZOS, 'v0', rand(N,1));
		kpt_vec = S.kptgrid(kpt,:);
		[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,kpt_vec);
		Hfun = @(x) h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,x,S,kpt_vec,spin);
        
        if ~(isreal(DL11) && isreal(DL22) && isreal(DL33)) || S.SOC_flag
            opts.isreal = false;
        end

        if(count == 1)
            % Upper bound estimator
            opts.maxit = 300; % WARNING: might need more accuracy
            if(S.ForceCount == 1)
                % For first relaxation step
                [upper_bound_guess_vecs(:,ks), bup(ks)] = (eigs(Hfun,N,1,'lr',opts));
                bup(ks) = real(bup(ks)) * 1.01;
            else
                % For subsequent relaxation steps
                opts.v0 = upper_bound_guess_vecs(:,ks);
                [upper_bound_guess_vecs(:,ks), bup(ks)] = eigs(Hfun,N,1,'lr',opts) ;
                bup(ks) = real(bup(ks)) * 1.01;
            end
            % Lower bound estimator
            if(S.ForceCount == 1)
                % For first relaxation step
                a0(ks) = real(eigs(Hfun,N,1,'sr',opts)) - 0.1;
            else
                % For subsequent relaxation steps use lowest eigenvalue
                a0(ks) = min(EigVal(:,ks));
            end
            % Always use this on every first SCF of relaxation
            %lambda_cutoff(ks) = 0.5 * (bup(ks) + a0(ks));
        else
            % For subsequent steps
            if count > S.rhoTrigger
                if S.chefsibound_flag == 1 || ((S.xc == 4) && (count == S.rhoTrigger + 1))
                    % for metaGGA: since 1st SCF of SCAN is GGA_PBE, it is
                    % necessary to make another Lanczos in 2nd SCF
                    % Upper bound
                    opts.tol = S.TOL_LANCZOS; % WARNING: might need more accuracy than the default
                    opts.v0 = upper_bound_guess_vecs(:,ks);
                    [upper_bound_guess_vecs(:,ks), bup(ks)] = eigs(Hfun,N,1,'lr',opts);
                    bup(ks) = real(bup(ks)) * 1.01;
                end
                % Lower bound
                a0(ks) = min(EigVal(:,ks));
            end
            % Set the filter cut off
            %S.lambda_f + log10(1e6-1)/S.bet + 0.5; 
            %lambda_cutoff(ks) = max(EigVal(:,ks)) + 0.10; 
        end

        if count == 1 && S.ForceCount == 1
            lambda_cutoff(ks) = 0.5 * (bup(ks) + a0(ks));
        else
            % Set the filter cut off
            %S.lambda_f + log10(1e6-1)/S.bet + 0.5; 
            lambda_cutoff(ks) = max(EigVal(:,ks)) + 0.10; 
        end

        % fprintf('filter cutoff = %f, lower bound = %f, upper bound = %f\n',lambda_cutoff(ks),a0(ks),bup(ks));
        % Chebyshev filtering
        psi(:,:,ks) = chebyshev_filter(psi(:,:,ks),S.npl,lambda_cutoff(ks),bup(ks),a0(ks),DL11,DL22,DL33,DG1,DG2,DG3,Veff,S,kpt_vec,spin);

        if S.StandardEigenFlag == 1
            psi(:,:,ks) = orth(psi(:,:,ks));
            Nev1 = size(psi(:,:,ks),2);    % WARNING: ORTH(psi) might change the size of psi, that's why we update Nev
            assert(Nev1 == S.Nev,'Number of states have changed within SCF');

            % Subspace Hamiltonian
            Hs = psi(:,:,ks)' * h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,psi(:,:,ks),S,kpt_vec,spin);

            % Solve subspace eigenproblem,
            if S.cell_typ < 3
                Hs = 0.5 * (Hs + Hs');
            end
            [Q, Q1] = eig(Hs);
            EigVal(:,ks) = real(diag(Q1)); % WARNING: Taking real part only!

            % subspace rotation
            psi(:,:,ks) = psi(:,:,ks) * Q;
        else
            % Subspace Hamiltonian
            Hs = psi(:,:,ks)' * h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,psi(:,:,ks),S,kpt_vec,spin);
            Ms = psi(:,:,ks)' * psi(:,:,ks);

            % Solve subspace eigenproblem
            if S.cell_typ < 3
                Hs = 0.5 * (Hs + Hs');
                Ms = 0.5 * (Ms + Ms');
            end

            [Q, Q1] = eig(Hs,Ms);
            EigVal(:,ks) = real(diag(Q1)); % WARNING: Taking real part only!
            psi(:,:,ks) = psi(:,:,ks) * Q;
        end
        
        % Normalize psi, s.t. integral(psi_new' * psi_new) = 1
        scfac = 1 ./ sqrt(sum(repmat(S.W,S.nspinor_eig,S.Nev) .* (psi(:,:,ks) .* conj(psi(:,:,ks))),1));
        % psi(:,:,ks) = psi(:,:,ks) * diag(scfac);
        psi(:,:,ks) = bsxfun(@times, psi(:,:,ks), scfac);
    end
else 
	% Before getting into parfor, set to use only one thread
	LASTN = maxNumCompThreads(1);
    parfor (ks = 1:S.tnkpt*S.nspin, S.num_worker_heuristic)
        if ks <= S.tnkpt
            kpt = ks;
            spin = 1;
        else
            kpt = ks - S.tnkpt;
            spin = 2;
        end
        % get Veff 
        if S.spin_typ == 2
            Veff = S.Veff;
        else
            Veff = S.Veff(:,spin);
        end

        rng('default'); % Initialize random number generator
		rng(ks+1);
		opts = struct('maxit', 100, 'tol', S.TOL_LANCZOS, 'v0', rand(N,1));
		kpt_vec = S.kptgrid(kpt,:);
		[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,kpt_vec);
		Hfun = @(x) h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,x,S,kpt_vec,spin);
        
        if ~(isreal(DL11) && isreal(DL22) && isreal(DL33)) || S.SOC_flag
            opts.isreal = false;
        end

        if(count == 1)
            % Upper bound estimator
            opts.maxit = 300; % WARNING: might need more accuracy
            if(S.ForceCount == 1)
                % For first relaxation step
                [upper_bound_guess_vecs(:,ks), bup(ks)] = (eigs(Hfun,N,1,'lr',opts));
                bup(ks) = real(bup(ks)) * 1.01;
            else
                % For subsequent relaxation steps
                opts.v0 = upper_bound_guess_vecs(:,ks);
                [upper_bound_guess_vecs(:,ks), bup(ks)] = eigs(Hfun,N,1,'lr',opts) ;
                bup(ks) = real(bup(ks)) * 1.01;
            end
            % Lower bound estimator
            if(S.ForceCount == 1)
                % For first relaxation step
                a0(ks) = real(eigs(Hfun,N,1,'sr',opts)) - 0.1;
            else
                % For subsequent relaxation steps use lowest eigenvalue
                a0(ks) = min(EigVal(:,ks));
            end
            % Always use this on every first SCF of relaxation
            %lambda_cutoff(ks) = 0.5 * (bup(ks) + a0(ks));
        else
            % For subsequent steps
            if count > S.rhoTrigger
                if S.chefsibound_flag == 1 || ((S.xc == 4) && (count == S.rhoTrigger + 1))
                    % for metaGGA: since 1st SCF of SCAN is GGA_PBE, it is
                    % necessary to make another Lanczos in 2nd SCF
                    % Upper bound
                    opts.tol = S.TOL_LANCZOS; % WARNING: might need more accuracy than the default
                    opts.v0 = upper_bound_guess_vecs(:,ks);
                    [upper_bound_guess_vecs(:,ks), bup(ks)] = eigs(Hfun,N,1,'lr',opts);
                    bup(ks) = real(bup(ks)) * 1.01;
                end
                % Lower bound
                a0(ks) = min(EigVal(:,ks));
            end
            % Set the filter cut off
            %S.lambda_f + log10(1e6-1)/S.bet + 0.5; 
            %lambda_cutoff(ks) = max(EigVal(:,ks)) + 0.10; 
        end

        if count == 1 && S.ForceCount == 1
            lambda_cutoff(ks) = 0.5 * (bup(ks) + a0(ks));
        else
            % Set the filter cut off
            %S.lambda_f + log10(1e6-1)/S.bet + 0.5; 
            lambda_cutoff(ks) = max(EigVal(:,ks)) + 0.10; 
        end

        % fprintf('filter cutoff = %f, lower bound = %f, upper bound = %f\n',lambda_cutoff(ks),a0(ks),bup(ks));
        % Chebyshev filtering
        psi(:,:,ks) = chebyshev_filter(psi(:,:,ks),S.npl,lambda_cutoff(ks),bup(ks),a0(ks),DL11,DL22,DL33,DG1,DG2,DG3,Veff,S,kpt_vec,spin);

        if S.StandardEigenFlag == 1
            psi(:,:,ks) = orth(psi(:,:,ks));
            Nev1 = size(psi(:,:,ks),2);    % WARNING: ORTH(psi) might change the size of psi, that's why we update Nev
            assert(Nev1 == S.Nev,'Number of states have changed within SCF');

            % Subspace Hamiltonian
            Hs = psi(:,:,ks)' * h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,psi(:,:,ks),S,kpt_vec,spin);

            % Solve subspace eigenproblem,
            if S.cell_typ < 3
                Hs = 0.5 * (Hs + Hs');
            end
            [Q, Q1] = eig(Hs);
            EigVal(:,ks) = real(diag(Q1)); % WARNING: Taking real part only!

            % subspace rotation
            psi(:,:,ks) = psi(:,:,ks) * Q;
        else
            % Subspace Hamiltonian
            Hs = psi(:,:,ks)' * h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,psi(:,:,ks),S,kpt_vec,spin);
            Ms = psi(:,:,ks)' * psi(:,:,ks);

            % Solve subspace eigenproblem
            if S.cell_typ < 3
                Hs = 0.5 * (Hs + Hs');
                Ms = 0.5 * (Ms + Ms');
            end

            [Q, Q1] = eig(Hs,Ms);
            EigVal(:,ks) = real(diag(Q1)); % WARNING: Taking real part only!
            psi(:,:,ks) = psi(:,:,ks) * Q;
        end
        
        % Normalize psi, s.t. integral(psi_new' * psi_new) = 1
        scfac = 1 ./ sqrt(sum(repmat(S.W,S.nspinor_eig,S.Nev) .* (psi(:,:,ks) .* conj(psi(:,:,ks))),1));
        % psi(:,:,ks) = psi(:,:,ks) * diag(scfac);
        psi(:,:,ks) = bsxfun(@times, psi(:,:,ks), scfac);
    end
	% After parfor, reset the number of threads as before
	maxNumCompThreads(LASTN);
end

% restore the storage of psi, EigVal
if S.spin_typ == 1
    psi = cat(1,psi(:,:,1:S.tnkpt),psi(:,:,1+S.tnkpt:2*S.tnkpt));    
end
end