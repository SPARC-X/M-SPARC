function stress_exx = evaluateHybridStress(S)
%**********************************************************************
%*                   Stress contribution from exact exchange          *
%**********************************************************************
stress_exx = zeros(3,3);
diag_term = 0;
if S.exxdivmethod == 1
    diag_term = S.Eex/4;
end

for spinor = 1:S.nspinor
    ndrange = (1+(spinor-1)*S.N:spinor*S.N); 
    nsshift = (spinor-1)*S.tnkpt*(S.spin_typ == 1);

    for k_ind = 1:S.tnkpt
        for q_ind = 1:S.tnkpthf
            % q_ind_rd is the index in reduced kptgrid
            q_ind_rd = S.kpthf_ind(q_ind,1);
            for i = 1:S.Nev
                for j = 1:S.Nev
                    if S.kpthf_ind(q_ind,2)
                        psiqi = S.psi(ndrange,i,q_ind_rd);
                    else
                        psiqi = conj(S.psi(ndrange,i,q_ind_rd));
                    end
                    psikj = S.psi(ndrange,j,k_ind);
                    rhs = conj(psiqi) .* psikj;

                    k = S.kptgrid(k_ind,:);
                    q = S.kptgridhf(q_ind,:);
                    k_shift = k - q;
                    phi1 = poissonSolve_FFT(S,rhs,k_shift,S.const_stress);
                    if S.exxdivmethod == 0
                        phi2 = poissonSolve_FFT(S,rhs,k_shift,S.const_stress_2);
                    end

                    Dphi_x = blochGradient(S,k_shift,1)*phi1;
                    Dphi_y = blochGradient(S,k_shift,2)*phi1;
                    Dphi_z = blochGradient(S,k_shift,3)*phi1;

                    Dcrho_x = conj(blochGradient(S,k_shift,1)*rhs);
                    Dcrho_y = conj(blochGradient(S,k_shift,2)*rhs);
                    Dcrho_z = conj(blochGradient(S,k_shift,3)*rhs);

                    stress_exx(1,1) = stress_exx(1,1) - S.wkpt(k_ind)*S.wkpthf(q_ind)*S.occ_outer(i,q_ind_rd+nsshift)*S.occ_outer(j,k_ind+nsshift)*real(sum(S.hyb_mixing.*Dcrho_x.*Dphi_x.*S.W));
                    stress_exx(2,2) = stress_exx(2,2) - S.wkpt(k_ind)*S.wkpthf(q_ind)*S.occ_outer(i,q_ind_rd+nsshift)*S.occ_outer(j,k_ind+nsshift)*real(sum(S.hyb_mixing.*Dcrho_y.*Dphi_y.*S.W));
                    stress_exx(3,3) = stress_exx(3,3) - S.wkpt(k_ind)*S.wkpthf(q_ind)*S.occ_outer(i,q_ind_rd+nsshift)*S.occ_outer(j,k_ind+nsshift)*real(sum(S.hyb_mixing.*Dcrho_z.*Dphi_z.*S.W));
                    stress_exx(1,2) = stress_exx(1,2) - S.wkpt(k_ind)*S.wkpthf(q_ind)*S.occ_outer(i,q_ind_rd+nsshift)*S.occ_outer(j,k_ind+nsshift)*real(sum(S.hyb_mixing.*Dcrho_x.*Dphi_y.*S.W));
                    stress_exx(1,3) = stress_exx(1,3) - S.wkpt(k_ind)*S.wkpthf(q_ind)*S.occ_outer(i,q_ind_rd+nsshift)*S.occ_outer(j,k_ind+nsshift)*real(sum(S.hyb_mixing.*Dcrho_x.*Dphi_z.*S.W));
                    stress_exx(2,3) = stress_exx(2,3) - S.wkpt(k_ind)*S.wkpthf(q_ind)*S.occ_outer(i,q_ind_rd+nsshift)*S.occ_outer(j,k_ind+nsshift)*real(sum(S.hyb_mixing.*Dcrho_y.*Dphi_z.*S.W));
                    if S.exxdivmethod == 0
                        diag_term = diag_term - S.wkpt(k_ind)*S.wkpthf(q_ind)*S.occ_outer(i,q_ind_rd+nsshift)*S.occ_outer(j,k_ind+nsshift)*real(sum(S.hyb_mixing.*conj(rhs).*phi2.*S.W));
                    end
                end
            end
        end
    end
end

stress_exx(2,1) = stress_exx(1,2);
stress_exx(3,1) = stress_exx(1,3);
stress_exx(3,2) = stress_exx(2,3);
stress_exx = stress_exx/2*S.occfac;    

% convert to cartesian coordinates
stress_exx = S.grad_T'*stress_exx*S.grad_T;

% compute final stress_exx
stress_exx = 2*stress_exx + (2*diag_term-2*S.Eex)*eye(3); 

end
