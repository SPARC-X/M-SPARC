function force_hub = atomicForce_hubbard(S)
% @brief    Calculates the atomic force for each atom with U correction.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param S      Struct containing relevant fields
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================

fprintf('\nStarting hubbard force calculation ... \n');

% Calculate local forces
tic_hubForces = tic;

% Initialize
force_hub = zeros(S.n_atm,3);

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

for kpt = 1 : S.tnkpt
    if (kpt(1) == 0 && kpt(2) == 0 && kpt(3) == 0)
        fac = 1.0;
    else
        fac = 1.0i;
    end
    kpt_vec = S.kptgrid(kpt,:);

    % Calculate gradient of psi
    Dpsi_x = zeros(S.nspinor *S.N,S.Nev);
    Dpsi_y = zeros(S.nspinor *S.N,S.Nev);
    Dpsi_z = zeros(S.nspinor *S.N,S.Nev);
    for spinor = 1:S.nspinor
        ndrange = (1+(spinor-1)*S.N:spinor*S.N);
        Dpsi_x(ndrange,:) = blochGradient(S,kpt_vec,1)*S.psi(ndrange,:,kpt);
        Dpsi_y(ndrange,:) = blochGradient(S,kpt_vec,2)*S.psi(ndrange,:,kpt);
        Dpsi_z(ndrange,:) = blochGradient(S,kpt_vec,3)*S.psi(ndrange,:,kpt);
    end

    for spinor = 1:S.nspinor
        % shift = (spinor-1)*S.N;         % for selecting each spinor, spinor=1 shift = 0, spinor=2, shift = S.N
        nsshift = (spinor-1)*S.tnkpt*(S.spin_typ == 1);
        ndrange = (1+(spinor-1)*S.N:spinor*S.N);

        % scalar relativistic part
        force_atm = zeros(S.n_atm,3);

        % psi
        psi = S.psi(ndrange,:,kpt);

        % Loop over all atoms having U correction
        for JJ_a = 1 : S.n_atm_U

            Orb = zeros(size(psi,1),S.AtomU(JJ_a).angnum);
            for img = 1:S.AtomU(JJ_a).n_image_rc
                Orb(S.AtomU(JJ_a).rcImage(img).rc_pos,:) = Orb(S.AtomU(JJ_a).rcImage(img).rc_pos,:)...
                    + S.AtomU(JJ_a).rcImage(img).orb_mat * (exp(-dot(kpt_vec,(S.AtomsU(JJ_a,:)-S.AtomU(JJ_a).rcImage(img).coordinates)*fac)));
            end

            wt_Orb = bsxfun(@times,Orb,S.W);

            rho_mn = S.AtomU(JJ_a).rho_mn(:,:,spinor);
            pre_fac = (0.5*eye(size(rho_mn)) - rho_mn); % Size: m x m

            integral_1_x = (Dpsi_x(ndrange,:))' * wt_Orb; % Nev x m
            integral_1_y = (Dpsi_y(ndrange,:))' * wt_Orb; % Nev x m
            integral_1_z = (Dpsi_z(ndrange,:))' * wt_Orb; % Nev x m
            integral_2 = wt_Orb'*psi; % m x Nev

            sum_1_x = integral_2 * (S.occ(:,kpt+nsshift) .* integral_1_x); % m x m
            sum_1_y = integral_2 * (S.occ(:,kpt+nsshift) .* integral_1_y); % m x m
            sum_1_z = integral_2 * (S.occ(:,kpt+nsshift) .* integral_1_z); % m x m

            tf_x = (S.AtomU(JJ_a).U) .* (pre_fac * 2 * real(sum_1_x)); % m x m
            tf_x = trace(tf_x); % Summing over m

            tf_y = (S.AtomU(JJ_a).U)' .* (pre_fac * 2 * real(sum_1_y)); % m x m
            tf_y = trace(tf_y); % Summing over m

            tf_z = (S.AtomU(JJ_a).U)' .* (pre_fac * 2 * real(sum_1_z)); % m x m
            tf_z = trace(tf_z); % Summing over m

            idx = corresponding_indices(JJ_a);
            force_atm(idx , :) = force_atm(idx , :) + [tf_x tf_y tf_z];
        end
        force_hub = force_hub - S.occfac*S.wkpt(kpt)*force_atm;
    end

end
fprintf('Hubbard force calculation: %.3f s\n', toc(tic_hubForces));
end