function pres_exx = evaluateHybridPressure(S)
pres_exx = 0;
if S.usefock > 0 && S.exxdivmethod ~= 0
    for spin = 1:S.nspin
        spin_shift = (spin-1)*S.tnkpt;
        for k_ind = 1:S.tnkpt
            for q_ind = 1:S.tnkpthf
                % q_ind_rd is the index in reduced kptgrid
                q_ind_rd = S.kpthf_ind(q_ind,1);
                for i = 1:S.Nev
                    for j = 1:S.Nev
                        if S.kpthf_ind(q_ind,2)
                            psiqi = S.psi(:,i,q_ind_rd+spin_shift);
                        else
                            psiqi = conj(S.psi(:,i,q_ind_rd+spin_shift));
                        end
                        psikj = S.psi(:,j,k_ind+spin_shift);
                        rhs = conj(psiqi) .* psikj;

                        k = S.kptgrid(k_ind,:);
                        q = S.kptgridhf(q_ind,:);
                        k_shift = k - q;
                        phi = poissonSolve_FFT(S,rhs,k_shift,S.const_press);
                        pres_exx = pres_exx - S.wkpt(k_ind)*S.wkpthf(q_ind)*S.occ_outer(i,q_ind_rd+spin_shift)*S.occ_outer(j,k_ind+spin_shift)*real(sum(S.hyb_mixing.*conj(rhs).*phi.*S.W));
                    end
                end
            end
        end
    end
    pres_exx = pres_exx/2*S.occfac;
    if S.exxdivmethod == 1
        pres_exx = pres_exx + 3/4*S.Eex;
    end
    pres_exx = 2*pres_exx - 2*S.Eex;
end
end