function AtomU = calculate_loc_orbitals(S)
% @brief    Calculates local orbital projectors for each atom with U correction.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param S      Struct containing relevant fields
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================

t1 = tic;
fprintf('\n Starting calculating local orbitals ...\n');

% Initialization
AtomU = repmat(struct([]),[1,S.n_atm_U]);

% Orbital occupation calculation
count_U_typ = 1;
count_U_typ_atms = 1;

% loop over all atoms with U correction
for JJ_a = 1:S.n_atm_U
    % Load info about isolated orbitals
    lmax = S.AtmU(count_U_typ).lmax;
    AtomU(JJ_a).lmax = lmax;
    AtomU(JJ_a).count_U_typ = count_U_typ;

    % Atom position of atom JJ_a
    x0 = S.AtomsU(JJ_a,1);
    y0 = S.AtomsU(JJ_a,2);
    z0 = S.AtomsU(JJ_a,3);

    % Store cutoff zone for local orbitals
    if (S.cell_typ == 1)
        rcx = S.AtmU(count_U_typ).rcU;
        rcy = S.AtmU(count_U_typ).rcU;
        rcz = S.AtmU(count_U_typ).rcU;
    else
        rcx = S.AtmU(count_U_typ).rbU_x;
        rcy = S.AtmU(count_U_typ).rbU_y;
        rcz = S.AtmU(count_U_typ).rbU_z;
    end

    % Calculate number of images in each direction based on b/c
    if S.BCx == 0
        n_image_xl = floor((S.AtomsU(JJ_a,1) + rcx)/S.L1);
        n_image_xr = floor((S.L1 - S.AtomsU(JJ_a,1)+rcx-S.dx)/S.L1);
    else
        n_image_xl = 0;
        n_image_xr = 0;
    end

    if S.BCy == 0
        n_image_yl = floor((S.AtomsU(JJ_a,2) + rcy)/S.L2);
        n_image_yr = floor((S.L2 - S.AtomsU(JJ_a,2)+rcy-S.dy)/S.L2);
    else
        n_image_yl = 0;
        n_image_yr = 0;
    end

    if S.BCz == 0
        n_image_zl = floor((S.AtomsU(JJ_a,3) + rcz)/S.L3);
        n_image_zr = floor((S.L3 - S.AtomsU(JJ_a,3)+rcz-S.dz)/S.L3);
    else
        n_image_zl = 0;
        n_image_zr = 0;
    end

    % Total no. of images of atom JJ_a (including atom JJ_a)
    n_image_total = (n_image_xl+n_image_xr+1) * (n_image_yl+n_image_yr+1) * (n_image_zl+n_image_zr+1);

    % Find the coordinates for all the images
    xx_img = (-n_image_xl : n_image_xr) * S.L1 + x0;
    yy_img = (-n_image_yl : n_image_yr) * S.L2 + y0;
    zz_img = (-n_image_zl : n_image_zr) * S.L3 + z0;
    [XX_IMG_3D,YY_IMG_3D,ZZ_IMG_3D] = ndgrid(xx_img,yy_img,zz_img);

    % Count for rc-images (images whoes rc-domain have influence in fund. domain)
    imgRcCount = 0;

    % Loop over all image(s) of atom JJ_a (including atom JJ_a)
    for count_image = 1:n_image_total
        % Atom position of the image
        x0_i = XX_IMG_3D(count_image);
        y0_i = YY_IMG_3D(count_image);
        z0_i = ZZ_IMG_3D(count_image);


        % Check if the rc-region of the image has influence in the fund. domain
        isImgRc_x = ((x0_i-S.xin)-S.L1+S.dx <= rcx) || (abs(x0_i-S.xin) <= rcx);
        isImgRc_y = ((y0_i-S.yin)-S.L2+S.dy <= rcy) || (abs(y0_i-S.yin) <= rcy);
        isImgRc_z = ((z0_i-S.zin)-S.L3+S.dz <= rcz) || (abs(z0_i-S.zin) <= rcz);
        isImgRc = isImgRc_x && isImgRc_y && isImgRc_z;
        if(isImgRc)
            % Starting and ending indices of orbital rc-region
            ii_rc_s = ceil( (x0_i - rcx)/S.dx) + 1;
            ii_rc_e = floor((x0_i + rcx)/S.dx) + 1;
            jj_rc_s = ceil( (y0_i - rcy)/S.dy) + 1;
            jj_rc_e = floor((y0_i + rcy)/S.dy) + 1;
            kk_rc_s = ceil( (z0_i - rcz)/S.dz) + 1;
            kk_rc_e = floor((z0_i + rcz)/S.dz) + 1;

            ii_rc_s = max(ii_rc_s,1);
            ii_rc_e = min(ii_rc_e,S.Nx);
            jj_rc_s = max(jj_rc_s,1);
            jj_rc_e = min(jj_rc_e,S.Ny);
            kk_rc_s = max(kk_rc_s,1);
            kk_rc_e = min(kk_rc_e,S.Nz);

            dd_nl = zeros(ii_rc_e-ii_rc_s+1, jj_rc_e-jj_rc_s+1, kk_rc_e-kk_rc_s+1);
            XX = dd_nl;
            YY = dd_nl;
            ZZ = dd_nl;
            rc_pos = dd_nl;
            rc_pos_ii = dd_nl;
            rc_pos_jj = dd_nl;
            rc_pos_kk = dd_nl;
            count3 = 1;
            for kk=kk_rc_s:kk_rc_e
                zz = S.zin + (kk - 1) * S.dz;
                count2 = 1;
                for jj=jj_rc_s:jj_rc_e
                    yy = S.yin + (jj - 1) * S.dy;
                    count1 = 1;
                    for ii=ii_rc_s:ii_rc_e
                        xx = S.xin + (ii - 1) * S.dx;
                        XX(count1,count2,count3) = xx;% - x0_i;
                        YY(count1,count2,count3) = yy;% - y0_i;
                        ZZ(count1,count2,count3) = zz;% - z0_i;
                        row_count = (kk-1)*S.Nx*S.Ny + (jj-1)*S.Nx + ii;
                        rc_pos(count1,count2,count3) = row_count;
                        rc_pos_ii(count1,count2,count3) = ii;
                        rc_pos_jj(count1,count2,count3) = jj;
                        rc_pos_kk(count1,count2,count3) = kk;
                        count1 = count1 + 1;
                    end
                    count2 = count2 + 1;
                end
                count3 = count3 + 1;
            end
            dd_nl = calculateDistance(XX,YY,ZZ,x0_i,y0_i,z0_i,S);
            dd_nl = reshape(dd_nl,[],1);
            rc_pos_ii = reshape(rc_pos_ii,[],1);
            rc_pos_jj = reshape(rc_pos_jj,[],1);
            rc_pos_kk = reshape(rc_pos_kk,[],1);
            XX = reshape(XX,[],1);
            YY = reshape(YY,[],1);
            ZZ = reshape(ZZ,[],1);
            rc_pos = reshape(rc_pos,[],1);

            isInBall = (dd_nl <= S.AtmU(count_U_typ).rcU);
            % Find the ball with radius rc
            rc_dd_nl = dd_nl(isInBall);
            rc_n_nodes = length(rc_dd_nl);
            if(rc_n_nodes~=0)
                imgRcCount = imgRcCount + 1;
                XX = XX(isInBall); YY = YY(isInBall); ZZ = ZZ(isInBall);
                pos_node_cart = coordinateTransformation(S,[XX,YY,ZZ],'noncart2cart_dis');
                pos_atm_cart = coordinateTransformation(S,[x0_i,y0_i,z0_i],'noncart2cart_dis');
                rc_XX = bsxfun(@minus,pos_node_cart(:,1),pos_atm_cart(:,1));
                rc_YY = bsxfun(@minus,pos_node_cart(:,2),pos_atm_cart(:,2));
                rc_ZZ = bsxfun(@minus,pos_node_cart(:,3),pos_atm_cart(:,3));

                AtomU(JJ_a).rcImage(imgRcCount).coordinates = [x0_i y0_i z0_i];
                AtomU(JJ_a).rcImage(imgRcCount).rc_pos = rc_pos(isInBall);
                AtomU(JJ_a).rcImage(imgRcCount).rc_pos_ii = rc_pos_ii(isInBall);
                AtomU(JJ_a).rcImage(imgRcCount).rc_pos_jj = rc_pos_jj(isInBall);
                AtomU(JJ_a).rcImage(imgRcCount).rc_pos_kk = rc_pos_kk(isInBall);
                AtomU(JJ_a).rcImage(imgRcCount).rc_n_nodes = rc_n_nodes;

                angnum = 0;
                count_l = 1; % to keep track of count of l
                for l = 0 : lmax % for all s, p, d, f states
                    if S.AtmU(count_U_typ).Uval(l+1) ~= 0 % Only for non-zero U value perform this
                        % Storing states: s - 2, p - 3, d - 5, f - 7
                        if l == 0
                            states_of_l = 1;
                        elseif l == 1
                            states_of_l = 3;
                        elseif l == 2
                            states_of_l = 5;
                        elseif l == 3
                            states_of_l = 7;
                        end
                        if S.AtmU(count_U_typ).hub_core == 1
                            l_strt = 1; % If want all l-th orbitals
                        else
                            l_strt = S.AtmU(count_U_typ).nmax_l(l+1); % If only want valence
                        end
                        for l_count = l_strt:S.AtmU(count_U_typ).nmax_l(l+1) 
                            for m = -l:l
                                angnum = angnum + 1;
                                AtomU(JJ_a).l_occ(angnum) = S.AtmU(count_U_typ).orb(l+1).count(l_count);
                                AtomU(JJ_a).l_states(angnum) = states_of_l;
                            end
                            count_l = count_l + 1;
                        end
                    end

                end

                AtomU(JJ_a).angnum = angnum;
                AtomU(JJ_a).rcImage(imgRcCount).orb_mat = zeros(rc_n_nodes,AtomU(JJ_a).angnum);
                count_orb = 1;
                for l = 0:lmax % for s, p, d and f states
                    if S.AtmU(count_U_typ).Uval(l+1) ~= 0 % Only for non-zero U value perform this
                        if S.AtmU(count_U_typ).hub_core == 1
                            l_strt = 1; % If want all l-th orbitals
                        else
                            l_strt = S.AtmU(count_U_typ).nmax_l(l+1); % If only want valence
                        end
                        for l_count = l_strt:S.AtmU(count_U_typ).nmax_l(l+1) % Comment out this for statement if you want only in the valence
                            orbital = interp1(S.AtmU(count_U_typ).r, S.AtmU(count_U_typ).orb(l+1).proj(:,l_count),rc_dd_nl,'spline');
                            for m = -l:l
                                Ylm = sphericalHarmonics(rc_XX,rc_YY,rc_ZZ,l,m,'real');
                                AtomU(JJ_a).rcImage(imgRcCount).orb_mat(:,count_orb) = AtomU(JJ_a).rcImage(imgRcCount).orb_mat(:,count_orb) + orbital .* Ylm;
                                count_orb = count_orb + 1;
                            end
                        end
                    end
                end
            else
                isImgRc = false;
            end
        end % end of if
        AtomU(JJ_a).imgRb(count_image).isImgRc = isImgRc;
    end % end of loop over images (including atom JJ_a)

    AtomU(JJ_a).U = zeros(AtomU(JJ_a).angnum,1);
    count_orb = 1;
    for l = 0:lmax % for all s, p, d and f states
        if S.AtmU(count_U_typ).Uval(l+1) ~= 0 % Only for non-zero U value perform this
            if S.AtmU(count_U_typ).hub_core == 1
                l_strt = 1; % If want all l-th orbitals
            else
                l_strt = S.AtmU(count_U_typ).nmax_l(l+1); % If only want valence
            end
            for l_count = l_strt:S.AtmU(count_U_typ).nmax_l(l+1) % Comment out this for statement if you want only in the valence
                for m=-l:l
                    AtomU(JJ_a).U(count_orb) = S.AtmU(count_U_typ).Uval(l+1);
                    count_orb = count_orb + 1;
                end
            end
        end
    end
    

    % number of rc-images: rc-region which has influence in the fund. domain
    AtomU(JJ_a).n_image_rc = imgRcCount;

    % Check if same type of atoms are over
    if count_U_typ_atms == S.AtmU(count_U_typ).n_atm_typ
        count_U_typ_atms = 1;
        count_U_typ = count_U_typ + 1;
    else
        count_U_typ_atms = count_U_typ_atms + 1;
    end
    AtomU(JJ_a).ImgRbcount = n_image_total;
end % end of loop over atoms with U correction

fprintf(' Done. (%f s)\n',toc(t1));
end