function stress = evaluateHubbardStress(S)
% @brief    Calculates the hubbard stress contribution to periodic cells.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param S      Struct containing relevant fields
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================

fprintf('\nStarting hubbard stress calculation ... \n');

% Calculate local forces
tic_hub = tic;

stress = zeros(3,3);
S_T = transpose(S.lat_uvec);

for spinor = 1 : S.nspinor
    for JJ_a = 1 : S.n_atm_U
        rho_mn = S.AtomU(JJ_a).rho_mn(:,:,spinor);
        pre_fac = (0.5*eye(size(rho_mn)) - rho_mn); % Size: m x m
        pre_fac = (S.AtomU(JJ_a).U) .* pre_fac;
        E_term = S.occfac*trace(pre_fac*rho_mn);

        stress(1,1) = stress(1,1) - E_term;
        stress(2,2) = stress(2,2) - E_term;
        stress(3,3) = stress(3,3) - E_term;
    end
end

for kpt = 1 : S.tnkpt
    fac = 1.0i;
    kpt_vec = S.kptgrid(kpt,:);

    Dpsi_x = zeros(S.nspinor *S.N,S.Nev);
    Dpsi_y = zeros(S.nspinor *S.N,S.Nev);
    Dpsi_z = zeros(S.nspinor *S.N,S.Nev);
    for spinor = 1:S.nspinor
        ndrange = (1+(spinor-1)*S.N:spinor*S.N); 
        Dpsi_x(ndrange,:) = blochGradient(S,kpt_vec,1)*S.psi(ndrange,:,kpt);
        Dpsi_y(ndrange,:) = blochGradient(S,kpt_vec,2)*S.psi(ndrange,:,kpt);
        Dpsi_z(ndrange,:) = blochGradient(S,kpt_vec,3)*S.psi(ndrange,:,kpt);
    end

    TDpsi_1 = S.grad_T(1,1)*Dpsi_x + S.grad_T(2,1)*Dpsi_y + S.grad_T(3,1)*Dpsi_z ;
    TDpsi_2 = S.grad_T(1,2)*Dpsi_x + S.grad_T(2,2)*Dpsi_y + S.grad_T(3,2)*Dpsi_z ;
    TDpsi_3 = S.grad_T(1,3)*Dpsi_x + S.grad_T(2,3)*Dpsi_y + S.grad_T(3,3)*Dpsi_z ;

    for spinor = 1 : S.nspinor
        nsshift = (spinor-1)*S.tnkpt*(S.spin_typ == 1);
        shift = (spinor-1)*S.N;         % for selecting each spinor, spinor=1 shift = 0, spinor=2, shift = S.N
        ndrange = (1+(spinor-1)*S.N:spinor*S.N);

        % psi
        psi = S.psi(ndrange,:,kpt);

        for JJ_a = 1 : S.n_atm_U
            rho_mn = S.AtomU(JJ_a).rho_mn(:,:,spinor);
            pre_fac = (0.5*eye(size(rho_mn)) - rho_mn); % Size: m x m
            pre_fac = (S.AtomU(JJ_a).U) .* pre_fac;

            Orb = zeros(size(psi,1),S.AtomU(JJ_a).angnum);

            integral_2_xx = zeros(S.AtomU(JJ_a).angnum,S.Nev);
            integral_2_xy = zeros(S.AtomU(JJ_a).angnum,S.Nev);
            integral_2_xz = zeros(S.AtomU(JJ_a).angnum,S.Nev);
            integral_2_yy = zeros(S.AtomU(JJ_a).angnum,S.Nev);
            integral_2_yz = zeros(S.AtomU(JJ_a).angnum,S.Nev);
            integral_2_zz = zeros(S.AtomU(JJ_a).angnum,S.Nev);

            for img = 1:S.AtomU(JJ_a).n_image_rc
                Orb(S.AtomU(JJ_a).rcImage(img).rc_pos,:) = Orb(S.AtomU(JJ_a).rcImage(img).rc_pos,:)...
                    + S.AtomU(JJ_a).rcImage(img).orb_mat * (exp(-dot(kpt_vec,(S.AtomsU(JJ_a,:)-S.AtomU(JJ_a).rcImage(img).coordinates)*fac)));
            end

            wt_Orb = bsxfun(@times,Orb,S.W);

            integral_1 = conj(wt_Orb'*psi); % m x Nev

            for img = 1:S.AtomU(JJ_a).n_image_rc
                phase_fac = (exp(dot(kpt_vec,(S.AtomsU(JJ_a,:)-S.AtomU(JJ_a).rcImage(img).coordinates)*fac)));
                OrbW = transpose(bsxfun(@times, conj(S.AtomU(JJ_a).rcImage(img).orb_mat), S.W(S.AtomU(JJ_a).rcImage(img).rc_pos)));
                xr =(S.AtomU(JJ_a).rcImage(img).rc_pos_ii-1)*S.dx - S.AtomU(JJ_a).rcImage(img).coordinates(1) ;
                yr =(S.AtomU(JJ_a).rcImage(img).rc_pos_jj-1)*S.dy - S.AtomU(JJ_a).rcImage(img).coordinates(2) ;
                zr =(S.AtomU(JJ_a).rcImage(img).rc_pos_kk-1)*S.dz - S.AtomU(JJ_a).rcImage(img).coordinates(3) ;
                x_1 = S_T(1,1)*xr + S_T(1,2)*yr + S_T(1,3)*zr;
                y_1 = S_T(2,1)*xr + S_T(2,2)*yr + S_T(2,3)*zr;
                z_1 = S_T(3,1)*xr + S_T(3,2)*yr + S_T(3,3)*zr;

                integral_2_xx = integral_2_xx + OrbW * ...
                    ((TDpsi_1(S.AtomU(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(x_1,1,S.Nev)) * phase_fac ;

                integral_2_xy = integral_2_xy + OrbW * ...
                    ((TDpsi_1(S.AtomU(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(y_1,1,S.Nev)) * phase_fac ;

                integral_2_xz = integral_2_xz + OrbW * ...
                    ((TDpsi_1(S.AtomU(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(z_1,1,S.Nev)) * phase_fac ;

                integral_2_yy = integral_2_yy + OrbW * ...
                    ((TDpsi_2(S.AtomU(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(y_1,1,S.Nev)) * phase_fac ;

                integral_2_yz = integral_2_yz + OrbW * ...
                    ((TDpsi_2(S.AtomU(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(z_1,1,S.Nev)) * phase_fac ;

                integral_2_zz = integral_2_zz + OrbW * ...
                    ((TDpsi_3(S.AtomU(JJ_a).rcImage(img).rc_pos+shift,:)).*repmat(z_1,1,S.Nev)) * phase_fac ;
            end

            tf_xx = 2 * real(integral_2_xx * (S.occ(:,kpt+nsshift) .* transpose(integral_1)));
            tf_xy = 2 * real(integral_2_xy * (S.occ(:,kpt+nsshift) .* transpose(integral_1)));
            tf_xz = 2 * real(integral_2_xz * (S.occ(:,kpt+nsshift) .* transpose(integral_1)));
            tf_yy = 2 * real(integral_2_yy * (S.occ(:,kpt+nsshift) .* transpose(integral_1)));
            tf_yz = 2 * real(integral_2_yz * (S.occ(:,kpt+nsshift) .* transpose(integral_1)));
            tf_zz = 2 * real(integral_2_zz * (S.occ(:,kpt+nsshift) .* transpose(integral_1)));

            % Update stress tensor
            stress(1,1) = stress(1,1) - S.occfac * S.wkpt(kpt) * trace(pre_fac * tf_xx);
            stress(1,2) = stress(1,2) - S.occfac * S.wkpt(kpt) * trace(pre_fac * tf_xy);
            stress(1,3) = stress(1,3) - S.occfac * S.wkpt(kpt) * trace(pre_fac * tf_xz);
            stress(2,2) = stress(2,2) - S.occfac * S.wkpt(kpt) * trace(pre_fac * tf_yy);
            stress(2,3) = stress(2,3) - S.occfac * S.wkpt(kpt) * trace(pre_fac * tf_yz);
            stress(3,3) = stress(3,3) - S.occfac * S.wkpt(kpt) * trace(pre_fac * tf_zz);
        end
    end
end

cell_measure = S.Jacb;
if S.BCx == 0
	cell_measure = cell_measure * S.L1;
end
if S.BCy == 0
	cell_measure = cell_measure * S.L2;
end
if S.BCz == 0
	cell_measure = cell_measure * S.L3;
end

stress(2,1) = stress(1,2);
stress(3,1) = stress(1,3);
stress(3,2) = stress(2,3);

stress = stress / cell_measure;
fprintf('Hubbard stress calculation: %.3f s\n', toc(tic_hub));
end