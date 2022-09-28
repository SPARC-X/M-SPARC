function Atom = calculate_nloc_projector(S)
% @brief	calculate_nloc_projector(S) calculates nonlocal projectors for
%           each atom.
%
% @param S  A struct that contains the relevant fields.
%
% @authors	Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%			Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%

t1 = tic;
fprintf('\n Starting calculating nonlocal projectors ...\n');

% initialization
Atom = repmat(struct([]),[1,S.n_atm]);

% Pseudocharge generation and self energy calculation
count_typ = 1;
count_typ_atms = 1;

for JJ_a = 1:S.n_atm % loop over all the atoms
	% Load pseudopotential file and isolated atom electron density
	lloc = S.Atm(count_typ).lloc;
	lmax = S.Atm(count_typ).lmax;
	Atom(JJ_a).lloc = lloc;
	Atom(JJ_a).lmax = lmax;
	Atom(JJ_a).psptyp = S.Atm(count_typ).psptyp;
    Atom(JJ_a).count_typ = count_typ;
	if S.Atm(count_typ).psptyp == 1
		Atom(JJ_a).nproj = S.Atm(count_typ).nproj;
        if S.Atm(count_typ).pspsoc == 1
            Atom(JJ_a).nprojso = S.Atm(count_typ).nprojso;
        end
	end
	
	% Atom position of atom JJ_a
	x0 = S.Atoms(JJ_a,1);
	y0 = S.Atoms(JJ_a,2);
	z0 = S.Atoms(JJ_a,3);
	
	% Note the S.dx, S.dy, S.dz terms are to ensure the image rc-region overlap w/ fund. domain
	if (S.cell_typ == 1)
		rcx = S.Atm(count_typ).rc;
		rcy = S.Atm(count_typ).rc;
		rcz = S.Atm(count_typ).rc;        
	else
		rcx = S.Atm(count_typ).rb_x;
		rcy = S.Atm(count_typ).rb_y;
		rcz = S.Atm(count_typ).rb_z;
	end
	
	if S.BCx == 0
		n_image_xl = floor((S.Atoms(JJ_a,1) + rcx)/S.L1);
		n_image_xr = floor((S.L1 - S.Atoms(JJ_a,1)+rcx-S.dx)/S.L1);
	else
		n_image_xl = 0;
		n_image_xr = 0;
	end
	
	if S.BCy == 0
		n_image_yl = floor((S.Atoms(JJ_a,2) + rcy)/S.L2);
		n_image_yr = floor((S.L2 - S.Atoms(JJ_a,2)+rcy-S.dy)/S.L2);
	else
		n_image_yl = 0;
		n_image_yr = 0;
	end
	
	if S.BCz == 0
		n_image_zl = floor((S.Atoms(JJ_a,3) + rcz)/S.L3);
		n_image_zr = floor((S.L3 - S.Atoms(JJ_a,3)+rcz-S.dz)/S.L3);
	else
		n_image_zl = 0;
		n_image_zr = 0;
	end
	
	% Total No. of images of atom JJ_a (including atom JJ_a)
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

		% Check if the rc-region is inside the domain in Dirichlet BC
		% direction
		% isInside = (S.BCx == 0 || (S.BCx == 1 && (ii_s>1) && (ii_e<S.Nx))) && ...
		% 	(S.BCy == 0 || (S.BCy == 1 && (jj_s>1) && (jj_e<S.Ny))) && ...
		% 	(S.BCz == 0 || (S.BCz == 1 && (kk_s>1) && (kk_e<S.Nz)));
		% this might be redundant since rb is larger than rc
		% assert(isInside,'ERROR: Atom too close to boundary for nonlocal projector calculation');

		%**********************************************************************
		%*                   Calculate nonlocal projector Chi                 *
		%**********************************************************************
		% Check if the rc-region of the image has influence in the fund. domain
		isImgRc_x = ((x0_i-S.xin)-S.L1+S.dx <= rcx) || (abs(x0_i-S.xin) <= rcx);
		isImgRc_y = ((y0_i-S.yin)-S.L2+S.dy <= rcy) || (abs(y0_i-S.yin) <= rcy);
		isImgRc_z = ((z0_i-S.zin)-S.L3+S.dz <= rcz) || (abs(z0_i-S.zin) <= rcz);
		isImgRc = isImgRc_x && isImgRc_y && isImgRc_z;
		if(isImgRc)
			% Starting and ending indices of the nonlocal rc-region
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
			
			% if ~((ii_rc_e < ii_rc_s) || (jj_rc_e < jj_rc_s) || (kk_rc_e < kk_rc_s))
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

			isInBall = (dd_nl <= S.Atm(count_typ).rc);
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
				if (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
					fac1 = (y0-y0_i)/S.L2;
					fac2 = (z0-z0_i)/S.L3;
					ROT = (S.RotM1^fac1) * (S.RotM2^fac2);
					P = transpose(ROT * [rc_XX,rc_YY,rc_ZZ]');
					rc_XX = P(:,1); rc_YY = P(:,2); rc_ZZ = P(:,3);
				end
				Atom(JJ_a).rcImage(imgRcCount).coordinates = [x0_i y0_i z0_i];
				Atom(JJ_a).rcImage(imgRcCount).rc_pos = rc_pos(isInBall);
				Atom(JJ_a).rcImage(imgRcCount).rc_pos_ii = rc_pos_ii(isInBall);
				Atom(JJ_a).rcImage(imgRcCount).rc_pos_jj = rc_pos_jj(isInBall);
				Atom(JJ_a).rcImage(imgRcCount).rc_pos_kk = rc_pos_kk(isInBall);
				Atom(JJ_a).rcImage(imgRcCount).rc_n_nodes = rc_n_nodes;

				%loop over l (from 0 to lmax, except Vloc component), then 
				%for each l, loop over m (which has 2*l+1 values)
				if S.Atm(count_typ).psptyp == 0
					Atom(JJ_a).angnum = (lmax^2 + 2*lmax + 1) - (2*lloc + 1);
					Atom(JJ_a).rcImage(imgRcCount).Chi_mat = zeros(rc_n_nodes, Atom(JJ_a).angnum);
					count_pp = 1;
					for l=0:lmax
						if l~=lloc
							% Evaluating UdV_Jl = U_Jl * (V_Jl - V_J)
							UdV_Jl = interp1(S.Atm(count_typ).r_grid_vloc, S.Atm(count_typ).UdV(:,l+1),rc_dd_nl,'spline');
							for m=-l:l
								Ylm = sphericalHarmonics(rc_XX,rc_YY,rc_ZZ,l,m,'real');
								Atom(JJ_a).rcImage(imgRcCount).Chi_mat(:,count_pp) = Atom(JJ_a).rcImage(imgRcCount).Chi_mat(:,count_pp) + UdV_Jl .* Ylm;
								count_pp = count_pp + 1;
							end
						end
					end
				elseif S.Atm(count_typ).psptyp == 1
					angnum = 0;
					for l=0:lmax
						if l~=lloc
							for pp = 1:S.Atm(count_typ).nproj(l+1)
								for m=-l:l
									angnum = angnum + 1;
								end
							end
						end
					end
					Atom(JJ_a).angnum = angnum;
					Atom(JJ_a).rcImage(imgRcCount).Chi_mat = zeros(rc_n_nodes, Atom(JJ_a).angnum);
					count_pp = 1;
					for l=0:lmax
						if l~=lloc
							for pp = 1:S.Atm(count_typ).nproj(l+1)
								UdV_Jl= interp1(S.Atm(count_typ).r_grid_vloc, S.Atm(count_typ).Pot(l+1).proj(:,pp),rc_dd_nl,'spline');
								for m=-l:l
									Ylm = sphericalHarmonics(rc_XX,rc_YY,rc_ZZ,l,m,'real');
									Atom(JJ_a).rcImage(imgRcCount).Chi_mat(:,count_pp) = Atom(JJ_a).rcImage(imgRcCount).Chi_mat(:,count_pp) + UdV_Jl .* Ylm;
									count_pp = count_pp + 1;
								end
							end
						end
					end
                end
                
                if S.Atm(count_typ).pspsoc == 1
                    % implementation for psp8 format
                    angnumso = 0;
					for l=1:lmax
						if l~=lloc
							for pp = 1:S.Atm(count_typ).nprojso(l)
								for m=-l:l
									angnumso = angnumso + 1;
								end
							end
						end
					end
					Atom(JJ_a).angnumso = angnumso;
					Atom(JJ_a).rcImage(imgRcCount).Chiso_mat = zeros(rc_n_nodes, Atom(JJ_a).angnumso);
					count_pp = 1;
					for l=1:lmax
						if l~=lloc
							for pp = 1:S.Atm(count_typ).nprojso(l)
								UdV_Jl= interp1(S.Atm(count_typ).r_grid_vloc, S.Atm(count_typ).Potso(l).proj(:,pp),rc_dd_nl,'spline');
								for m=-l:l
									Ylm = sphericalHarmonics(rc_XX,rc_YY,rc_ZZ,l,m,'complex');
									Atom(JJ_a).rcImage(imgRcCount).Chiso_mat(:,count_pp) = Atom(JJ_a).rcImage(imgRcCount).Chiso_mat(:,count_pp) + UdV_Jl .* Ylm;
									count_pp = count_pp + 1;
								end
							end
						end
					end
                end

			else
				isImgRc = false;
			end
		end % end if
		Atom(JJ_a).imgRb(count_image).isImgRc = isImgRc;
	end % end of loop over images (including atom JJ_a)
	
	Atom(JJ_a).gamma_Jl = zeros(Atom(JJ_a).angnum,1);
	count_pp = 1;
	if Atom(JJ_a).psptyp == 0
		for l=0:lmax
			if l~=lloc
				for m=-l:l
					Atom(JJ_a).gamma_Jl(count_pp) = gamma_Jl(l+1);
					count_pp = count_pp + 1;
				end
			end
		end
	elseif  Atom(JJ_a).psptyp == 1
		for l=0:lmax
			if l~=lloc
				for pp = 1:S.Atm(count_typ).nproj(l+1)
					for m=-l:l
						Atom(JJ_a).gamma_Jl(count_pp) = S.Atm(count_typ).Pot(l+1).gamma_Jl(pp);
						count_pp = count_pp + 1;
					end
				end
			end
		end
    end
	
    if S.Atm(count_typ).pspsoc == 1
        % gamma for first term 0.5*m*gamma_Jln
        Atom(JJ_a).term1_index_so = zeros(Atom(JJ_a).angnumso,1);
        Atom(JJ_a).term1_gammaso_Jl = zeros(Atom(JJ_a).angnumso,1);
        indx_so = 0;
        count_pp = 0;
        for l=1:lmax
            if l~=lloc
                for pp = 1:S.Atm(count_typ).nprojso(l)
                    gammaso_Jl = S.Atm(count_typ).Potso(l).gamma_Jl(pp);
                    for m=-l:l
                        indx_so = indx_so + 1;
                        % only store the projectors with nonzero scaled factor
                        if m ~= 0
                            count_pp = count_pp + 1;
                            Atom(JJ_a).term1_index_so(count_pp) = indx_so;
                            Atom(JJ_a).term1_gammaso_Jl(count_pp,1) = 0.5*m*gammaso_Jl;
                        end
                    end
                end
            end
        end
        Atom(JJ_a).ncol_term1 = count_pp;
        
        % gamma and indexing for second term
        Atom(JJ_a).term2_index_so = zeros(Atom(JJ_a).angnumso,1);
        Atom(JJ_a).term2_gammaso_Jl = zeros(Atom(JJ_a).angnumso,1);
        indx_so = 0;
        count_pp = 0;
        for l=1:lmax
            if l~=lloc
                for pp = 1:S.Atm(count_typ).nprojso(l)
                    gammaso_Jl = S.Atm(count_typ).Potso(l).gamma_Jl(pp);
                    for m=-l:l
                        indx_so = indx_so + 1;
                        % only store the projectors with nonzero scaled factor
                        if m < l
                            count_pp = count_pp + 1;
                            Atom(JJ_a).term2_index_so(count_pp) = indx_so;
                            Atom(JJ_a).term2_gammaso_Jl(count_pp) = sqrt(l*(l+1)-m*(m+1))*gammaso_Jl*0.5;
                        end
                    end
                end
            end
        end
        Atom(JJ_a).ncol_term2 = count_pp;
    end
    
	% number of rc-images: rc-region of which has influence in the fund. domain
	Atom(JJ_a).n_image_rc = imgRcCount;

	% Check if same type of atoms are over
	if count_typ_atms == S.Atm(count_typ).n_atm_typ
		count_typ_atms = 1;
		count_typ = count_typ + 1;
	else
		count_typ_atms = count_typ_atms + 1;
	end
	Atom(JJ_a).ImgRbcount = n_image_total;
end % end of loop over atoms

fprintf(' Done. (%f s)\n',toc(t1));

end
