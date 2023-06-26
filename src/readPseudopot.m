function S = readPseudopot(S, ityp, psdfname, element)
% @brief    READPSEUDOPOT reads the pseudopotential file (psp8 format).
%
% @param ityp       Element type index.
% @param psdfname   The pseudopotential filename, with suffix.
% @param element    Element type.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%

%filename = sprintf('%s/%s',S.inputfile_path, psdfname);
filename = fullfile(S.inputfile_path, psdfname);

fid = fopen(filename,'r') ;

assert(fid~=-1,'Error: Cannot open pseudopotential file %s',filename);

textscan(fid,'%s',1,'delimiter','\n') ;
fscanf(fid,'%s',1);
Z = fscanf(fid,'%f',1);
textscan(fid,'%s',1,'delimiter','\n');
A = textscan(fid,'%f %f %f %f %f %f') ;

lmax = A{1,3};
lloc = A{1,4}; 
% lloc = 4; % forcing nonlocal component to 4!
mmax = A{1,5};

textscan(fid,'%s',1,'delimiter','\n') ;
Ax = textscan(fid,'%f %f %f') ;
textscan(fid,'%s',1,'delimiter','\n') ;
fchrg = Ax{2};

A = textscan(fid,'%f %f %f %f %f');
nproj = ones(lmax+1,1);
for i = 0:lmax
	nproj(i+1) = A{1,i+1};
end
textscan(fid,'%s',1,'delimiter','\n') ;

extension_switch = fscanf(fid,'%f',1);
textscan(fid,'%s',1,'delimiter','\n') ;

nprojso = zeros(lmax,1);
pspsoc = 0; % indicating if the psp file including spin-orbit coupling
if extension_switch == 2 || extension_switch == 3
    fprintf("This psp8 includes spin-orbit coupling.\n");
    pspsoc = 1;
    A = textscan(fid,'%f %f %f %f %f');
    for i = 1:lmax
        nprojso(i) = A{1,i};
    end
    textscan(fid,'%s',1,'delimiter','\n');
end

Pot = repmat(struct([]), lmax+1, 1);
Potso = repmat(struct([]), lmax, 1);

l_read = fscanf(fid,'%f',1);
for l = 0:lmax
	if l ~= lloc
		A = fscanf(fid,'%f',nproj(l+1)) ;
		Pot(l+1).gamma_Jl = A(:); 
		sz = [2+nproj(l+1),mmax];
		A = fscanf(fid,'%g',sz) ;
		r = A(2,:)' ;
		Pot(l+1).proj = A(3:end,:)';
		Pot(l+1).proj(2:end,:) = Pot(l+1).proj(2:end,:)./repmat(r(2:end),1,nproj(l+1));
		Pot(l+1).proj(1,:) = Pot(l+1).proj(2,:);
	else
		textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line  
		A = fscanf(fid,'%g',[3,mmax]);
		r = A(2,:)' ;
		Vloc = A(3,:)';
	end
	l_read = fscanf(fid,'%f',1);
end

if lloc > lmax || l_read == 4
	A = fscanf(fid,'%g',[3,mmax]) ;
	r = A(2,:)';
	Vloc = A(3,:)';
else
	% move back file pointer 4 columns
	fseek(fid, -4, 'cof');
end

% read spin-orbit projectors
if pspsoc == 1
    for l = 1:lmax
        fscanf(fid,'%f',1);
        A = fscanf(fid,'%f',nprojso(l)) ;
        Potso(l).gamma_Jl = A(:); 
        sz = [2+nprojso(l),mmax];
        A = fscanf(fid,'%g',sz) ;
        r = A(2,:)' ;
        Potso(l).proj = A(3:end,:)';
        Potso(l).proj(2:end,:) = Potso(l).proj(2:end,:)./repmat(r(2:end),1,nprojso(l));
        Potso(l).proj(1,:) = Potso(l).proj(2,:);
    end
end

% read core density
if fchrg > 0
	uu = zeros(mmax,1);
	Atilde = fscanf(fid,'%d %g %g %g %g %g %g',[7,mmax]) ;
	uu(1:end,1) = Atilde(3,:)/(4*pi);
	rho_Tilde = uu;
	rTilde = Atilde(2,:)';
    S.NLCC_flag = 1;
else
	rTilde = r;
	rho_Tilde = zeros(length(r),1);
end

uu = zeros(mmax,1);
asize = [5,mmax];
A = fscanf(fid,'%d %g %g %g %g',asize);
uu(1:end,1) = A(3,:)/(4*pi);
rho_isolated_guess = uu;

frewind(fid);
fscanf(fid,'%s',3);
rc = 0 ;
for l = 0:lmax
	r_core_read = fscanf(fid,'%g',1);
    rc_max = r_core_read;
    if l ~= lloc
        % check if r_core is large enough s.t. |proj| < 1E-8
        r_indx = find(r < r_core_read,1,'last');
        for i = 1:size(Pot(l+1).proj,2)
            try
                rc_temp = r(r_indx + find(abs(Pot(l+1).proj(r_indx+1:end,i)) < 1E-8,1) - 1);
            catch
                rc_temp = r(end);
            end
            if rc_temp > rc_max
                rc_max = rc_temp;
            end
        end
        fprintf("atom type %d, l = %d, r_core read %.5f, change to rmax where |UdV| < 1E-8, %.5f.\n", ityp, l, r_core_read, rc_max);
    end
    if rc_max > rc 
		rc = rc_max;
    end
end

r_grid_vloc = r;
r_grid_rho = r;

S.Atm(ityp).Z = Z;
S.Atm(ityp).Vloc = Vloc;
S.Atm(ityp).r_grid_vloc = r_grid_vloc;
S.Atm(ityp).rc = rc; 
S.Atm(ityp).Pot = Pot;
S.Atm(ityp).lmax = lmax;
S.Atm(ityp).lloc = lloc;
S.Atm(ityp).nproj = nproj;
S.Atm(ityp).r_grid_rho = r_grid_rho;
S.Atm(ityp).rho_isolated_guess = rho_isolated_guess;
S.Atm(ityp).rho_Tilde = rho_Tilde;
S.Atm(ityp).r_grid_rho_Tilde = rTilde;
S.Atm(ityp).pspsoc = pspsoc;
S.Atm(ityp).Potso = Potso;
S.Atm(ityp).nprojso = nprojso;

fclose(fid);
end
