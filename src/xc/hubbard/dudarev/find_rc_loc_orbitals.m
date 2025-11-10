function S = find_rc_loc_orbitals(S)
% @brief    Sets the cutoff radius for the local orbitals for atoms with U
%           correction.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param S      Struct containing relevant fields
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================

n_Uatm_typ = 1;
rc = 0;
for ityp = 1:S.n_typ
    if S.U_typ_flags(ityp) == 1
        % Find rc for where projector < tolerance
        lmax = S.AtmU(n_Uatm_typ).lmax;
        for l = 0:lmax
            if S.AtmU(n_Uatm_typ).Uval(l+1) == 0
                continue;
            end
            rc_max = S.Atm(ityp).rc_max(l+1);
            r_core = S.Atm(ityp).rc_max(l+1);
            r = S.AtmU(n_Uatm_typ).r;
            r_len = length(r);
            for l_count = 1:S.AtmU(n_Uatm_typ).nmax_l(l+1)
                for JJ = r_len : -1 : 1
                    if abs(S.AtmU(n_Uatm_typ).orb(l+1).proj(JJ,l_count)) > 5e-3
                        rc_temp = r(JJ);
                        break;
                    end
                end
                if rc_temp > rc_max
                    rc_max = rc_temp;
                end
            end
            if rc_max > rc
                rc = rc_max;
            end
        end

        n_Uatm_typ = n_Uatm_typ + 1;
    end
    if n_Uatm_typ > S.n_typ_U
        break;
    end
end

rc = ceil(rc);
fprintf("The local orbital cutoff is %f.\n", rc);

n_Uatm_typ = 1;
for ityp = 1:S.n_typ
    if S.U_typ_flags(ityp) == 1 
        S.AtmU(n_Uatm_typ).rcU = rc;

        % rb_U calc for non-ortho cell
        if S.cell_typ == 2
            % rb_x = S.AtmU(n_Uatm_typ).rcU;
        	% rb_y = S.AtmU(n_Uatm_typ).rcU;
        	% rb_z = S.AtmU(n_Uatm_typ).rcU;
            % rb_x = ceil(rb_x/S.dx-1e-12)*S.dx;
            % rb_y = ceil(rb_y/S.dy-1e-12)*S.dy;
            % rb_z = ceil(rb_z/S.dz-1e-12)*S.dz;
            LV = S.lat_vec;
            LV_inv = (pinv(LV))';
            norm_vec = [norm(LV_inv(1,:)) norm(LV_inv(2,:)) norm(LV_inv(3,:))];
            norm_inv_vec = 1./norm_vec;
            num_bins = ceil(S.AtmU(n_Uatm_typ).rcU./norm_inv_vec);
            rb_x = num_bins(1)*norm(LV(1,:));
            rb_y = num_bins(2)*norm(LV(2,:));
            rb_z = num_bins(3)*norm(LV(3,:));

            S.AtmU(n_Uatm_typ).rbU_x = rb_x;
            S.AtmU(n_Uatm_typ).rbU_y = rb_y;
            S.AtmU(n_Uatm_typ).rbU_z = rb_z;
        end

        n_Uatm_typ = n_Uatm_typ + 1;
    end
    if n_Uatm_typ > S.n_typ_U
        break;
    end
end

end