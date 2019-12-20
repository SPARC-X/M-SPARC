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
lloc = A{1,4}; lloc = 4; % forcing nonlocal component to 4!
mmax = A{1,5};

textscan(fid,'%s',2,'delimiter','\n') ;
A = textscan(fid,'%f %f %f %f %f');
nproj = ones(lmax+1,1);
for i = 0:lmax
	nproj(i+1) = A{1,i+1};
end
textscan(fid,'%s',2,'delimiter','\n') ;

Pot = repmat(struct([]), lmax+1, 1);
l = 0;
while l<=lmax
	fscanf(fid,'%f',1) ;
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
		A = fscanf(fid,'%g',[3,mmax]);
		r = A(2,:)' ;
		Vloc = A(3,:)';
	end
	l = l+1;
end

if lloc > lmax
	fscanf(fid,'%f',1) ;
	A = fscanf(fid,'%g',[3,mmax]) ;
	r = A(2,:)';
	Vloc = A(3,:)';
end

uu = zeros(mmax,1);
size = [5,mmax];
A = fscanf(fid,'%d %g %g %g %g',size);
uu(1:end,1) = A(3,:)/(4*pi);
rho_isolated_guess = uu;

frewind(fid);
fscanf(fid,'%s',3);
rc = 0 ;
for i=1:lmax+1
	a = fscanf(fid,'%g',1) ;
	if a > rc 
		rc = a ;
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

end
