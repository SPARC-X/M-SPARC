function  Hx = h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,X,S,kptvec,spin)
% @brief   Calculates Hamiltonian vector product, where Vnl is the nonlocal
%          pseudopotential to be determined using the info stored in S.
%
% @authors  Abhiraj Sharma <asharma424@gatech.edu>
%           Qimen Xu <qimenxu@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% 
% @param DLii        Discrete laplacian component in 1D along ith direction
% @param DGi         Discrete gradient component in 1D along ith direction
% @param Veff        Effective potential N x 1 vector
% @param X           N x Ns matrix of eigenfunctions of Hamiltonian
% @param kptvec      k-point vector for the current Block diagonalized problem
% @param Hx          Hamiltonian times vector  
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%========================================================================================

Hx = zeros(size(X));
if (kptvec(1) == 0 && kptvec(2) == 0 && kptvec(3) == 0) && (S.SOC_flag == 0)
    fac = 1.0;
else
    fac = 1.0i;
end

for spinor = 1:S.nspinor_eig
    ndrange = (1+(spinor-1)*S.N:spinor*S.N);
    ndrange_opo = (1+(2-spinor)*S.N:(3-spinor)*S.N);
    sigma = (-1)^(spinor-1);
    
    % apply Veff
    if S.spin_typ == 2
        Hx(ndrange,:) = Veff(:,spinor).*X(ndrange,:) + (Veff(:,3)+sigma*1i*Veff(:,4)).*X(ndrange_opo,:); 
    else
        Hx(ndrange,:) = Veff.*X(ndrange,:);
    end

    shift = (spinor-1)*S.N;         % for selecting each spinor, spinor=1 shift = 0, spinor=2, shift = S.N
    shift2 = (2-spinor)*S.N;        % for selecting the other spin channel, spinor=1 shift2 = S.N, spinor=2,shift2=0  

    Hx(ndrange,:) = Hx(ndrange,:) - 0.5*(lapVec(DL11,DL22,DL33,DG1,DG2,DG3,X(ndrange,:),S));
    
    % apply scalar relativistic part
    for J = 1:S.n_atm
        Chisc_X_mult = zeros(S.Atom(J).angnum,size(X,2));
        for img = 1:S.Atom(J).n_image_rc
            Chisc_X_mult = Chisc_X_mult + ( bsxfun(@times, S.Atom(J).rcImage(img).Chi_mat, S.W(S.Atom(J).rcImage(img).rc_pos)) )' * X(S.Atom(J).rcImage(img).rc_pos+shift,:) * ...
                        (exp(dot(kptvec,(S.Atoms(J,:)-S.Atom(J).rcImage(img).coordinates)*fac)));
        end
        Chisc_X_mult = bsxfun(@times,Chisc_X_mult, S.Atom(J).gamma_Jl);
        for img = 1:S.Atom(J).n_image_rc
            Chisc = S.Atom(J).rcImage(img).Chi_mat * (exp(-dot(kptvec,(S.Atoms(J,:)-S.Atom(J).rcImage(img).coordinates)*fac)));
            Hx(S.Atom(J).rcImage(img).rc_pos+shift,:) = Hx(S.Atom(J).rcImage(img).rc_pos+shift,:) + Chisc * Chisc_X_mult;
        end
    end

    % apply spin orbit part1
    for J = 1:S.n_atm
        if S.Atm(S.Atom(J).count_typ).pspsoc == 0
            continue;
        end
        ncol_term1 = S.Atom(J).ncol_term1;
        soindx = S.Atom(J).term1_index_so(1:ncol_term1);
        Chiso_X_mult = zeros(ncol_term1,size(X,2));
        for img = 1:S.Atom(J).n_image_rc
            Chiso_X_mult = Chiso_X_mult + ( bsxfun(@times, S.Atom(J).rcImage(img).Chiso_mat(:,soindx), S.W(S.Atom(J).rcImage(img).rc_pos)) )' * X(S.Atom(J).rcImage(img).rc_pos+shift,:) * ...
                        (exp(dot(kptvec,(S.Atoms(J,:)-S.Atom(J).rcImage(img).coordinates)*fac)));
        end

        Chiso_X_mult = bsxfun(@times,Chiso_X_mult, sigma*S.Atom(J).term1_gammaso_Jl(1:ncol_term1));
        for img = 1:S.Atom(J).n_image_rc
            Chiso = S.Atom(J).rcImage(img).Chiso_mat(:,soindx) * (exp(-dot(kptvec,(S.Atoms(J,:)-S.Atom(J).rcImage(img).coordinates)*fac)));
            Hx(S.Atom(J).rcImage(img).rc_pos+shift,:) = Hx(S.Atom(J).rcImage(img).rc_pos+shift,:) + Chiso * Chiso_X_mult;
        end
    end

    % apply spin orbit part2
    for J = 1:S.n_atm
        if S.Atm(S.Atom(J).count_typ).pspsoc == 0
            continue;
        end
        ncol_term2 = S.Atom(J).ncol_term2;
        Chiso_Jlmp1n_X_mult = zeros(ncol_term2,size(X,2));
        if spinor == 1
            soindx1 = S.Atom(J).term2_index_so(1:ncol_term2)+1;
            soindx2 = S.Atom(J).term2_index_so(1:ncol_term2);
        else 
            soindx1 = S.Atom(J).term2_index_so(1:ncol_term2);
            soindx2 = S.Atom(J).term2_index_so(1:ncol_term2)+1;
        end
        for img = 1:S.Atom(J).n_image_rc
            Chiso_Jlmp1n_X_mult = Chiso_Jlmp1n_X_mult + ( bsxfun(@times, S.Atom(J).rcImage(img).Chiso_mat(:,soindx1), S.W(S.Atom(J).rcImage(img).rc_pos)) )' * X(S.Atom(J).rcImage(img).rc_pos+shift2,:) * ...
                        (exp(dot(kptvec,(S.Atoms(J,:)-S.Atom(J).rcImage(img).coordinates)*fac)));
        end

        Chiso_Jlmp1n_X_mult = bsxfun(@times,Chiso_Jlmp1n_X_mult, S.Atom(J).term2_gammaso_Jl(1:ncol_term2));
        for img = 1:S.Atom(J).n_image_rc
            Chiso_Jlmn = S.Atom(J).rcImage(img).Chiso_mat(:,soindx2) * (exp(-dot(kptvec,(S.Atoms(J,:)-S.Atom(J).rcImage(img).coordinates)*fac)));
            Hx(S.Atom(J).rcImage(img).rc_pos+shift,:) = Hx(S.Atom(J).rcImage(img).rc_pos+shift,:) + Chiso_Jlmn * Chiso_Jlmp1n_X_mult;
        end
    end

    % hybrid
    if S.usefock > 1
        Hx(ndrange,:) = evaluateExactExchangePotential(S,X(ndrange,:),Hx(ndrange,:),kptvec,spin);
    end

    % mgga
    if (S.xc == 4) && (S.countPotential > 0) % metaGGA, set a flag to seperate it from the 1st PBE SCF computation
        Hx(ndrange,:) = evaluateMggaPotential(S,X(ndrange,:),Hx(ndrange,:),kptvec,spin);
    end
end
end
