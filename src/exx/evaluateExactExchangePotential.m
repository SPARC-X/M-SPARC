function Vexx = evaluateExactExchangePotential(S,X,kptvec,spin)
if S.ACEFlag == 0
    spin_shift = (spin-1)*S.tnkpt;
    Vexx = zeros(S.N,size(X,2));
    V_guess = rand(S.N,1);
    for i = 1:size(X,2)
        for j = 1:S.Nev
            for q_ind = 1:S.tnkpthf
                % q_ind_rd is the index in reduced kptgrid
                q_ind_rd = S.kpthf_ind(q_ind,1);
                if S.kpthf_ind(q_ind,2)
                    psiq = S.psi_outer(:,j,q_ind_rd+spin_shift);
                else
                    psiq = conj(S.psi_outer(:,j,q_ind_rd +spin_shift));
                end

                rhs = conj(psiq).*X(:,i);
                if S.exxmethod == 0             % solving in fourier space
                    q = S.kptgridhf(q_ind,:);
                    k_shift = kptvec - q;
                    V_ji = poissonSolve_FFT(S,rhs,k_shift,S.const_by_alpha);
                else                            % solving in real space
                    f = poisson_RHS(S,rhs);
                    [V_ji, flag] = pcg(-S.Lap_std,-f,1e-8,1000,S.LapPreconL,S.LapPreconU,V_guess);
                    assert(flag==0);
                    V_guess = V_ji;
                end

                Vexx(:,i) = Vexx(:,i) - S.wkpthf(q_ind)*S.occ_outer(j,q_ind_rd+spin_shift)*V_ji.*psiq;
            end
        end
    end
else 
    if S.isgamma == 1
        col = 1+(spin-1)*S.Ns_occ(1):S.Ns_occ(1)+(spin-1)*S.Ns_occ(2);
        Xi_times_psi = (transpose(S.Xi(:,col))*X)*S.dV;
        Vexx = -S.Xi(:,col)*Xi_times_psi;
    else
        col = 1+(spin-1)*S.Ns_occ(1):S.Ns_occ(1)+(spin-1)*S.Ns_occ(2);
        k_ind = find(ismembertol(S.kptgrid,kptvec,1e-8,'ByRows',true))+0;
        Xi_times_psi = S.Xi(:,col,k_ind)'*X*(S.dV);
        Vexx = - S.Xi(:,col,k_ind)*Xi_times_psi;
    end
end
    
end
