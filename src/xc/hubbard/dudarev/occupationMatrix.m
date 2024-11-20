function S = occupationMatrix(S)
% @brief    Calculates the occupation matrix for each atom.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================

for J = 1 : S.n_atm_U
    S.AtomU(J).rho_mn = zeros(S.AtomU(J).angnum,S.AtomU(J).angnum,S.nspinor);
end

% First SCF step guess
if S.useHubbard == 1
    if S.spin_typ == 1
        % Initialize an array to hold corresponding indices of U atoms
        corresponding_indices = zeros(S.n_atm_U, 1);

        % match tolerance
        tolerance = 1e-6;

        for i = 1:S.n_atm_U
            % Find the index of the row in S.Atoms that matches the current row in S.AtomsU
            idx = find(all(abs(S.Atoms - S.AtomsU(i, :)) < tolerance, 2));

            % If a match is found, store the index; otherwise, store -1
            if ~isempty(idx)
                corresponding_indices(i) = idx;
            else
                corresponding_indices(i) = -1;  % No match found
            end
        end
    end
    

    for J = 1 : S.n_atm_U
        isMagnetic = 0;

        if S.spin_typ == 1
            atomMag = S.atomMag(corresponding_indices(J));
            if atomMag > 0
                isMagnetic = 1;
                maj_sp = 1; % up spin filled first
                min_sp = 2;
            elseif atomMag < 0
                isMagnetic = 1;
                maj_sp = 2; % dw spin filled first
                min_sp = 1;
            end
        end

        angnum = S.AtomU(J).angnum;
        
        if isMagnetic == 1 % magnetic
            for ang = 1 : angnum
                l_states = S.AtomU(J).l_states(ang); % Total number of mag states per s (1) or p (3) or d (5) or f (7) orbitals
                l_occ = S.AtomU(J).l_occ(ang); % Total occupation of that "l" state
                if l_occ > l_states
                    S.AtomU(J).rho_mn(ang, ang, maj_sp) = 1.0;
                    S.AtomU(J).rho_mn(ang, ang, min_sp) = (l_occ - l_states)/l_states;
                else
                    S.AtomU(J).rho_mn(ang, ang, maj_sp) = (l_occ - l_states)/l_states;
                end
            end
        else % non magnetic
            for spinor = 1 : S.nspinor
                for ang = 1 : angnum
                    S.AtomU(J).rho_mn(ang, ang, spinor) = S.AtomU(J).l_occ(ang)/S.AtomU(J).l_states(ang)/2;
                end
            end
        end


    end
    return;
end

for spinor = 1:S.nspinor_eig
    shift = (spinor-1)*S.N;
    for J = 1:S.n_atm_U
        for spinor_psi = 1 : S.nspinor
            rho_mn = zeros(S.AtomU(J).angnum,S.AtomU(J).angnum);
            ndrange = (1+(spinor_psi-1)*S.N:spinor_psi*S.N);
            % shift_psi = (spinor_psi-1)*S.N;
            nsshift = (spinor_psi - 1)*S.tnkpt*(S.spin_typ == 1);
            for kpt = 1 : S.tnkpt

                kptvec = S.kptgrid(kpt,:);
                if (kptvec(1) == 0 && kptvec(2) == 0 && kptvec(3) == 0) && (S.SOC_flag == 0)
                    fac = 1.0;
                else
                    fac = 1.0i;
                end

                psi = S.psi(ndrange,:,kpt);

                Orb = zeros(size(psi,1),S.AtomU(J).angnum);
                for img = 1:S.AtomU(J).n_image_rc
                    Orb(S.AtomU(J).rcImage(img).rc_pos+shift,:) = Orb(S.AtomU(J).rcImage(img).rc_pos+shift,:)...
                        + S.AtomU(J).rcImage(img).orb_mat * (exp(-dot(kptvec,(S.AtomsU(J,:)-S.AtomU(J).rcImage(img).coordinates)*fac)));
                end

                wt_psi = bsxfun(@times,psi,S.W);
                wt_Orb = bsxfun(@times,Orb,S.W);

                Orb_psi_mult = Orb'*wt_psi; % Size: m x Ns
                Psi_orb_mult = psi'*wt_Orb; % Size: Ns x m
                occ_Psi_orb_mult = S.occ(:,kpt+nsshift).*(Psi_orb_mult); % Size: Ns x m
                rho_mn = rho_mn + S.wkpt(kpt)*(Orb_psi_mult)*occ_Psi_orb_mult; % Size: m x m
            end
            rho_mn = 0.5*(rho_mn' + rho_mn); % Make rho hermitean
            S.AtomU(J).rho_mn(:,:,spinor_psi) = real(rho_mn);
        end
    end
end

end