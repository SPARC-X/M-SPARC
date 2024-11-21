function [] = print_Occ_mat(S)
% @brief    Prints the occupation matrix for each atom.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================

% Print
E_Hub = 0;
for J = 1 : S.n_atm_U
    fprintf(2,' ---------Atom %d---------\n',J);
    for spinor = 1 : S.nspinor
        rho_mn = S.AtomU(J).rho_mn(:,:,spinor);
        % Print
        if S.nspin > 1
            fprintf('SPIN %d\n',spinor);
        end

        % Print eigenvalues, total local occupation of the atom, eigenvectors
        fprintf('Trace[rho_mn]: %6.6f\n',S.occfac*trace(rho_mn));
        [evec,eval] = eig(rho_mn);
        fprintf('Eigenvalues:\n')
        eval_p = diag(eval);
        for r = 1 : length(eval_p)
            fprintf('%6.3f  ',eval_p(r));
        end
        fprintf('\n');

        fprintf('Eigenvectors:\n')
        [rows,cols] = size(evec);
        for r = 1 : rows
            for c = 1 : cols
                fprintf('%6.3f  ',evec(r,c));
            end
            fprintf('\n');
        end

        % Print rho_mn
        fprintf('Occupation matrix:\n')
        for r = 1 : size(rho_mn,1)
            for c = 1 : size(rho_mn,2)
                fprintf('%6.3f  ',rho_mn(r,c));
            end
            fprintf('\n');
        end
    end
end
end