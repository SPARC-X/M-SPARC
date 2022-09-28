function  Hnlx = h_nonlocal_vector_mult(DL11,DL22,DL33,DG1,DG2,DG3,Veff,X,S,kptvec,spin)
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
% @param Hnlx        Hamiltonian times vector  
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%========================================================================================

if S.nspinor == 1
    Hnlx = -0.5*(lapVec(DL11,DL22,DL33,DG1,DG2,DG3,X,S)) + bsxfun(@times,Veff,X);

    % Vnl * X
    if (kptvec(1) == 0 && kptvec(2) == 0 && kptvec(3) == 0)
        fac = 1.0;
    else
        fac = 1.0i;
    end

    for J = 1:S.n_atm
        Chi_X_mult = zeros(S.Atom(J).angnum,size(X,2));
        for img = 1:S.Atom(J).n_image_rc
            % Chi_X_mult = Chi_X_mult + (S.Atom(J).rcImage(img).Chi_mat .* S.W(S.Atom(J).rcImage(img).rc_pos))' * X(S.Atom(J).rcImage(img).rc_pos,:)* ...
            Chi_X_mult = Chi_X_mult + ( bsxfun(@times, S.Atom(J).rcImage(img).Chi_mat, S.W(S.Atom(J).rcImage(img).rc_pos)) )' * X(S.Atom(J).rcImage(img).rc_pos,:) * ...
                        (exp(dot(kptvec,(S.Atoms(J,:)-S.Atom(J).rcImage(img).coordinates)*fac)));
        end
        % Chi_X_mult = Chi_X_mult .* S.Atom(J).gamma_Jl;
        Chi_X_mult = bsxfun(@times,Chi_X_mult, S.Atom(J).gamma_Jl);
        for img = 1:S.Atom(J).n_image_rc
            Chi = S.Atom(J).rcImage(img).Chi_mat * (exp(-dot(kptvec,(S.Atoms(J,:)-S.Atom(J).rcImage(img).coordinates)*fac)));
            Hnlx(S.Atom(J).rcImage(img).rc_pos,:) = Hnlx(S.Atom(J).rcImage(img).rc_pos,:) + Chi * Chi_X_mult;
        end
    end
    
    % hybrid 
    if S.usefock > 1
        Vexx = evaluateExactExchangePotential(S,X,kptvec,spin);
        Hnlx = Hnlx + S.hyb_mixing*Vexx;
    end
    
    % mgga
    if (S.xc == 4) && (S.countPotential > 0) % metaGGA, set a flag to seperate it from the 1st PBE SCF computation
        if S.nspin == 1
            VxcScan3 = S.VxcScan3;
        else % spin polarization mGSGA
            VxcScan3 = S.VxcScan3(:, spin); 
        end
        nCol = size(X, 2);
        if S.cell_typ == 2 % unorthogonal cell
            lapc_T = [S.lapc_T(1,1), S.lapc_T(2,1), S.lapc_T(3,1);
                S.lapc_T(2,1), S.lapc_T(2,2), S.lapc_T(3,2);
                S.lapc_T(3,1), S.lapc_T(3,2), S.lapc_T(3,3)];
            v3grad1 = blochGradient(S,kptvec,1) *X; 
            v3grad2 = blochGradient(S,kptvec,2) *X; 
            v3grad3 = blochGradient(S,kptvec,3) *X; 

            v3gradpsiTheKpt = [v3grad1(:), v3grad2(:), v3grad3(:)];
            v3gradpsiTheKptMLapT = v3gradpsiTheKpt*lapc_T;
            v3gradpsiTheKptMLapT_1 = VxcScan3 .* reshape(v3gradpsiTheKptMLapT(:, 1), S.N, nCol);
            v3gradpsiTheKptMLapT_2 = VxcScan3 .* reshape(v3gradpsiTheKptMLapT(:, 2), S.N, nCol);
            v3gradpsiTheKptMLapT_3 = VxcScan3 .* reshape(v3gradpsiTheKptMLapT(:, 3), S.N, nCol);
            Hnlx = Hnlx - 0.5*(blochGradient(S,kptvec,1)*v3gradpsiTheKptMLapT_1 + blochGradient(S,kptvec,2)*v3gradpsiTheKptMLapT_2 + blochGradient(S,kptvec,3)*v3gradpsiTheKptMLapT_3);
        else % orthogonal cell
            v3grad1 = VxcScan3 .* (blochGradient(S,kptvec,1) *X);
            v3grad2 = VxcScan3 .* (blochGradient(S,kptvec,2) *X);
            v3grad3 = VxcScan3 .* (blochGradient(S,kptvec,3) *X);
            Hnlx = Hnlx - 0.5*(blochGradient(S,kptvec,1)*v3grad1 + blochGradient(S,kptvec,2)*v3grad2 + blochGradient(S,kptvec,3)*v3grad3);
        end
    end
    
elseif S.nspinor == 2
    
    assert(size(Veff,1) == S.N);
    assert(size(X,1) == 2*S.N);
    Hnlx = [Veff; Veff].*X;
    fac = 1.0i;
    
    for spinor = 1:S.nspinor
        ndrange = (1+(spinor-1)*S.N:spinor*S.N);
        sigma = (-1)^(spinor-1);
        shift = (spinor-1)*S.N;         % for selecting each spinor, spinor=1 shift = 0, spinor=2, shift = S.N
        shift2 = (2-spinor)*S.N;        % for selecting the other spin channel, spinor=1 shift2 = S.N, spinor=2,shift2=0  

        Hnlx(ndrange,:) = Hnlx(ndrange,:) - 0.5*(lapVec(DL11,DL22,DL33,DG1,DG2,DG3,X(ndrange,:),S));
        
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
                Hnlx(S.Atom(J).rcImage(img).rc_pos+shift,:) = Hnlx(S.Atom(J).rcImage(img).rc_pos+shift,:) + Chisc * Chisc_X_mult;
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
                Hnlx(S.Atom(J).rcImage(img).rc_pos+shift,:) = Hnlx(S.Atom(J).rcImage(img).rc_pos+shift,:) + Chiso * Chiso_X_mult;
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
                Hnlx(S.Atom(J).rcImage(img).rc_pos+shift,:) = Hnlx(S.Atom(J).rcImage(img).rc_pos+shift,:) + Chiso_Jlmn * Chiso_Jlmp1n_X_mult;
            end
        end

    end    
end



end
