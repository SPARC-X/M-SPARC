function S = ace_operator(S)
if S.exxmethod == 1
    V_guess = rand(S.N,1);
end

S.Ns_occ = zeros(1,2);

if S.isgamma == 1
    
    for spin = 1:S.nspin
        S.Ns_occ(spin) = min(sum(S.occ_outer(:,spin)>1e-6)+S.EXXACEVal_state,S.Nev);
    end
    S.Xi = zeros(S.N,sum(S.Ns_occ));    % For storage of W and Xi
    
    for spin = 1:S.nspin
        % spin_shift = (spin-1)*S.tnkpt;
        spin_xi_shift = (spin-1)*S.Ns_occ(1);
        Ns = S.Ns_occ(spin);
        rhs = zeros(S.N,Ns);
        for i = 1:Ns
            rhs(:,i:Ns) = bsxfun(@times,S.psi_outer(:,i:Ns,spin),S.psi_outer(:,i,spin));
            V_i = zeros(S.N,Ns);
            for j = i:Ns
                if (S.occ_outer(i,spin) + S.occ_outer(j,spin) > 1e-4)
                    if S.exxmethod == 0             % solving in fourier space
                        V_i(:,j) = poissonSolve_FFT(S,rhs(:,j),[0,0,0],S.const_by_alpha);
                    else                            % solving in real space
                        f = poisson_RHS(S,rhs(:,j));
                        [V_i(:,j), flag] = pcg(-S.Lap_std,-f,1e-8,1000,S.LapPreconL,S.LapPreconU,V_guess);
                        assert(flag==0);
                        V_guess = V_i(:,j);
                    end
                end
            end
            S.Xi(:,(i+1:Ns)+spin_xi_shift) = S.Xi(:,(i+1:Ns)+spin_xi_shift) - S.occ_outer(i,spin)*bsxfun(@times,S.psi_outer(:,i,spin),V_i(:,(i+1:Ns)));
            S.Xi(:,i+spin_xi_shift) = S.Xi(:,i+spin_xi_shift) - bsxfun(@times,S.psi_outer(:,(i:Ns),spin),V_i(:,(i:Ns))) * S.occ_outer((i:Ns),spin);
        end
        
        col = 1+(spin-1)*S.Ns_occ(1):S.Ns_occ(1)+(spin-1)*S.Ns_occ(2);
        M = (transpose(S.psi_outer(:,1:Ns,spin))*S.Xi(:,col))*S.dV;
        L = chol(-M); 
        S.Xi(:,col) = S.Xi(:,col) * inv(L); % Do it efficiently
    end
else
    for spin = 1:S.nspin
        col = (1:S.tnkpt)+(spin-1)*S.tnkpt;
        S.Ns_occ(spin) = min(max(sum(S.occ_outer(:,col) > 1e-6))+S.EXXACEVal_state,S.Nev);
    end
    
    S.Xi = zeros(S.N,sum(S.Ns_occ),S.tnkpt);    % For storage of W and Xi
    
    for spin = 1:S.nspin
        Ns = S.Ns_occ(spin);
        spin_shift = (spin-1)*S.tnkpt;
        xi_shift = (spin-1)*S.Ns_occ(1);
        for k_ind = 1:S.tnkpt
            for q_ind = 1:S.tnkpthf
                % q_ind_rd is the index in reduced kptgrid
                q_ind_rd = S.kpthf_ind(q_ind,1);

                k = S.kptgrid(k_ind,:);
                q = S.kptgridhf(q_ind,:);
                if S.kpthf_ind(q_ind,2)
                    psi_q_set = S.psi_outer(:,1:Ns,q_ind_rd+spin_shift);
                else
                    psi_q_set = conj(S.psi_outer(:,1:Ns,q_ind_rd+spin_shift));
                end

                for i = 1:Ns
                    psi_k = S.psi_outer(:,i,k_ind+spin_shift);
                    rhs = conj(psi_q_set) .* psi_k;
                    k_shift = k - q;
                    V_i = zeros(S.N,Ns);
                    for j = 1:Ns
                        if S.occ_outer(j,q_ind_rd+spin_shift) > 1e-6
                            if S.exxmethod == 0             % solving in fourier space
                                V_i(:,j) = poissonSolve_FFT(S,rhs(:,j),k_shift,S.const_by_alpha);
                            else                            % solving in real space
                                f = poisson_RHS(S,rhs(:,j));
                                [V_i(:,j), flag] = pcg(-S.Lap_std,-f,1e-8,1000,S.LapPreconL,S.LapPreconU,V_guess);
                                assert(flag==0);
                                V_guess = V_i(:,j);
                            end
                        end
                    end

                    S.Xi(:,i+xi_shift,k_ind) = S.Xi(:,i+xi_shift,k_ind) - S.wkpthf(q_ind)*(psi_q_set.*V_i)*S.occ_outer(1:Ns,q_ind_rd+spin_shift);
                end
            end
            
            col = 1+(spin-1)*S.Ns_occ(1):S.Ns_occ(1)+(spin-1)*S.Ns_occ(2);
            
            M = S.psi_outer(:,1:Ns,k_ind+spin_shift)'*S.Xi(:,col,k_ind)*S.dV;
            % to ensure M is hermitian
            M = 0.5*(M+M');
            L = chol(-M);
            S.Xi(:,col,k_ind) = S.Xi(:,col,k_ind) * inv(L); % Do it efficiently
        end
    end
end
end

