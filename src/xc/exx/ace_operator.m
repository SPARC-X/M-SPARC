function S = ace_operator(S)
if S.exxmethod == 1
    V_guess = rand(S.N,1);
end

S.Ns_occ = max(sum(S.occ_outer>1e-6));
S.Ns_occ = min(S.Ns_occ+S.EXXACEVal_state,S.Nev);
Ns = S.Ns_occ;

if S.isgamma == 1
    S.Xi = zeros(S.nspinor*S.N,Ns);    % For storage of W and Xi
    
    for spinor = 1:S.nspinor
        ndrange = (1+(spinor-1)*S.N:spinor*S.N); 
        nsshift = (spinor-1)*(S.spin_typ == 1);

        rhs = zeros(S.N,Ns);
        for i = 1:Ns
            rhs(:,i:Ns) = bsxfun(@times,S.psi_outer(ndrange,i:Ns),S.psi_outer(ndrange,i));
            V_i = zeros(S.N,Ns);
            for j = i:Ns
                if (S.occ_outer(i,1+nsshift) + S.occ_outer(j,1+nsshift) > 1e-4)
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
            S.Xi(ndrange,(i+1:Ns)) = S.Xi(ndrange,(i+1:Ns)) - S.occ_outer(i,1+nsshift)*bsxfun(@times,S.psi_outer(ndrange,i),V_i(:,(i+1:Ns)));
            S.Xi(ndrange,i) = S.Xi(ndrange,i) - bsxfun(@times,S.psi_outer(ndrange,(i:Ns)),V_i(:,(i:Ns))) * S.occ_outer((i:Ns),1+nsshift);
        end
        
        M = (transpose(S.psi_outer(ndrange,1:Ns))*S.Xi(ndrange,:))*S.dV;
        L = chol(-M); 
        S.Xi(ndrange,:) = S.Xi(ndrange,:) / L; % Do it efficiently
    end
else

    S.Xi = zeros(S.N*S.nspinor,Ns,S.tnkpt);    % For storage of W and Xi
    
    for spinor = 1:S.nspinor
        ndrange = (1+(spinor-1)*S.N:spinor*S.N); 
        nsshift = (spinor-1)*S.tnkpt*(S.spin_typ == 1);
        
        for k_ind = 1:S.tnkpt
            for q_ind = 1:S.tnkpthf
                % q_ind_rd is the index in reduced kptgrid
                q_ind_rd = S.kpthf_ind(q_ind,1);

                k = S.kptgrid(k_ind,:);
                q = S.kptgridhf(q_ind,:);
                if S.kpthf_ind(q_ind,2)
                    psi_q_set = S.psi_outer(ndrange,1:Ns,q_ind_rd);
                else
                    psi_q_set = conj(S.psi_outer(ndrange,1:Ns,q_ind_rd));
                end

                for i = 1:Ns
                    psi_k = S.psi_outer(ndrange,i,k_ind);
                    rhs = conj(psi_q_set) .* psi_k;
                    k_shift = k - q;
                    V_i = zeros(S.N,Ns);
                    for j = 1:Ns
                        if S.occ_outer(j,q_ind_rd+nsshift) > 1e-6
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

                    S.Xi(ndrange,i,k_ind) = S.Xi(ndrange,i,k_ind) - S.wkpthf(q_ind)*(psi_q_set.*V_i)*S.occ_outer((1:Ns),q_ind_rd+nsshift);
                end
            end
            
            
            M = S.psi_outer(ndrange,1:Ns,k_ind)'*S.Xi(ndrange,:,k_ind)*S.dV;
            % to ensure M is hermitian
            M = 0.5*(M+M');
            L = chol(-M);
            S.Xi(ndrange,:,k_ind) = S.Xi(ndrange,:,k_ind) / L; % Do it efficiently
        end
    end
end
end

