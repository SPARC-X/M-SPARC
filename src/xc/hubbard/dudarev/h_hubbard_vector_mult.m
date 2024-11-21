function Hx = h_hubbard_vector_mult(Hx, X, S, kptvec, spin)
% @brief    Calculates the Hamiltonian vector product, where V_U is the
%           hubbard potential to be determined using info stored in S.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param Hx         Hamiltonian times vector
% @param X          N x Ns matrix of eigenfunctions of Hamiltonian
% @param kptvec     k-point vector for the current Block diagonalized problem
% @param spin       Spin channel
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================

if (kptvec(1) == 0 && kptvec(2) == 0 && kptvec(3) == 0) && (S.SOC_flag == 0)
    fac = 1.0;
else
    fac = 1.0i;
end

for spinor = 1:S.nspinor_eig
    ndrange = (1+(spinor-1)*S.N:spinor*S.N);

    shift = (spinor-1)*S.N;         % for selecting each spinor, spinor=1 shift = 0, spinor=2, shift = S.N

    % apply the Hubbard potential
    for J = 1:S.n_atm_U
        Orb = zeros(size(X,1),S.AtomU(J).angnum);
        for img = 1:S.AtomU(J).n_image_rc
            Orb(S.AtomU(J).rcImage(img).rc_pos+shift,:) = Orb(S.AtomU(J).rcImage(img).rc_pos+shift,:)...
                + S.AtomU(J).rcImage(img).orb_mat * (exp(-dot(kptvec,(S.AtomsU(J,:)-S.AtomU(J).rcImage(img).coordinates)*fac)));
        end

        wt_X = bsxfun(@times,X(ndrange,:),S.W);

        Orb_X_mult = Orb'*wt_X; % Size: m x Ns

        % Calculate rho_mn
        rho_mn = S.AtomU(J).rho_mn(:,:,spin);

        % prefactor - matrix way
        pre_fac = (0.5*eye(size(rho_mn)) - rho_mn); % Size: m x m
        pre_fac = bsxfun(@times, pre_fac, S.AtomU(J).U);

        Orb_X_mult = pre_fac*Orb_X_mult;
        Hx(ndrange,:) = Hx(ndrange,:) + Orb * Orb_X_mult;

    end

end
end