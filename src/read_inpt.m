function S = read_inpt(S, filename) 
% @brief    READ_INPT reads the filename.inpt file.
%
% @param fname  The input filename, without suffix.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics % @2016-2019 (c) Georgia Institute of Technology. Mechanics Group, Georgia Tech
%

fprintf(' Reading .inpt file ...\n');

S.filename = filename;

[inputfile_path,~,~] = fileparts(filename);

if isempty(inputfile_path)
	S.inputfile_path = '.';
else
	S.inputfile_path = inputfile_path;
end

fid1=fopen(strcat(filename,'.inpt'),'r');
if (fid1 == -1) 
	error('\nCannot open file "%s.inpt"\n',filename);
end

% print out the .inpt file on screen for debug
fprintf('\n\n<INPT>\n');
fprintf('# $ cat %s.inpt',filename);
type(strcat(filename,'.inpt')); 
fprintf('<\\INPT>\n\n');

Flag_smear_typ = 0;
Flag_Temp = 0;
Flag_elecT = 0;
Flag_accuracy = 0;
Flag_ionT = 0;
Flag_eqT = 0;
%Flag_ionT_end = 0;
Flag_cell = 0;
Flag_latvec_scale = 0;
Flag_kptshift = 0;

while(~feof(fid1))
	C_inpt = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
	str = char(C_inpt{:}); % for R2016b and later, can use string()
	% skip commented lines starting by '#'
	%if (isempty(str) || str(1) == '#' || strcmp(str,'undefined'))
	if (isempty(str) || str(1) == '#')
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		%fprintf('skipping current line!\n');
		continue;
	end
	% check variable name and assign value
	if (strcmp(str,'NP_SPIN_PARAL:'))
		msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line           
	elseif (strcmp(str,'NP_KPOINT_PARAL:'))
		msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line           
	elseif (strcmp(str,'NP_BAND_PARAL:'))
		msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line           
	elseif (strcmp(str,'NP_DOMAIN_PARAL:'))
		msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line           
	elseif (strcmp(str,'NP_DOMAIN_PHI_PARAL:'))
		msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line           
	elseif (strcmp(str,'CELL:'))
		C_param = textscan(fid1,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
        Flag_cell = 1;
		S.L1 = C_param{1};
		S.L2 = C_param{2};
		S.L3 = C_param{3};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line           
    elseif (strcmp(str,'LATVEC_SCALE:'))
		C_param = textscan(fid1,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
        Flag_latvec_scale = 1;
        S.Flag_latvec_scale = 1;
		S.latvec_scale_x = C_param{1};
		S.latvec_scale_y = C_param{2};
		S.latvec_scale_z = C_param{3};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line       
	elseif (strcmp(str,'TWIST_ANGLE:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.alph = C_param{1}; % in radian/Bohr
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'LATVEC:'))
		C_param = textscan(fid1,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		S.lat_vec(1,:) = cell2mat(C_param);
		C_param = textscan(fid1,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		S.lat_vec(2,:) = cell2mat(C_param);
		C_param = textscan(fid1,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		S.lat_vec(3,:) = cell2mat(C_param);
	elseif (strcmp(str,'BOUNDARY_CONDITION:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.BC = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		fprintf('WARNING: "BOUNDARY_CONDITION" is obsolete, use "BC" instead!\n');
	elseif (strcmp(str,'BC:'))
		C_param = textscan(fid1,'%s %s %s',1,'delimiter',' ','MultipleDelimsAsOne',1); % read smearing type
		bcx = char(C_param{1});
		bcy = char(C_param{2});
		bcz = char(C_param{3});
		if bcx == 'p' || bcx == 'P'
			S.BCx = 0;
		elseif bcx == 'd' || bcx == 'D'
			S.BCx = 1;
		end

		if bcy == 'p' || bcy == 'P'
			S.BCy = 0;
		elseif bcy == 'd' || bcy == 'D'
			S.BCy = 1;
		elseif bcy == 'c' || bcy == 'C'
			S.BCy = 0;
			S.cell_typ = 3;
		end

		if bcz == 'p' || bcz == 'P'
			S.BCz = 0;
		elseif bcz == 'd' || bcz == 'D'
			S.BCz = 1;
		elseif bcz == 'h' || bcz == 'H'
			S.BCz = 0;
			if(S.cell_typ == 3)
				S.cell_typ = 5;
			else
				S.cell_typ = 4;
			end
		end
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'SPIN_TYP:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.spin_typ = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'FD_ORDER:'))
		%fscanf(fid1,'%f',S.order);
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.FDn = C_param{1};
		S.FDn = S.FDn / 2; % FDn is half the order
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'FD_GRID:'))
		C_param = textscan(fid1,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.Nx = C_param{1};
		S.Ny = C_param{2};
		S.Nz = C_param{3};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'ECUT:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.ecut = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MESH_SPACING:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.mesh_spacing = C_param{1}; % mesh spacing
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'KPOINT_GRID:'))
		C_param = textscan(fid1,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.nkpt = cell2mat(C_param);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'KPOINT_SHIFT:'))	
		C_param = textscan(fid1,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);	
		S.kptshift = cell2mat(C_param);	
        Flag_kptshift = 1;
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line    
	elseif (strcmp(str,'ELEC_TEMP_TYPE:'))
		Flag_smear_typ = Flag_smear_typ + 1;
		C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1); % read smearing type
		temp = char(C_param{:}); 
		if (strcmp(temp,'fd') || strcmp(temp,'fermi-dirac'))
			% currently only fermi-dirac is implemented in M-SPARC
			%fprintf('Fermi-Dirac smearing\n');
			S.elec_T_type = 0;
		elseif (strcmp(temp,'gaussian'))
			% currently only fermi-dirac is implemented in M-SPARC
			%error('Gaussian smearing is not implemented in M-SPARC!\n');
			S.elec_T_type = 1;
		else
			error(['\nCannot recognize electronic temperature (smearing) type: "%s"\n',...
				  'Available options: "fd" (or "fermi-dirac"), "gaussian" \n'], temp);
			%fprintf('Available options: "fd" (or "fermi-dirac"), "gaussian" \n');
		end
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'SMEARING:'))
		Flag_Temp = Flag_Temp + 1;
		Flag_elecT = Flag_elecT + 1;
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		smearing = C_param{1};
		S.bet = 1 / smearing;
		S.Temp = smearing / S.kB;
		%str = 'undefined';    % initialize str
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'BETA:'))
		Flag_Temp = Flag_Temp + 1;
		Flag_elecT = Flag_elecT + 1;
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.bet = C_param{1};
		S.Temp = 1 / (S.kB * S.bet);
		str = 'undefined';    % initialize str
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'ELEC_TEMP:'))
		Flag_Temp = Flag_Temp + 1;
		Flag_elecT = Flag_elecT + 1;
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.Temp = C_param{1};
		S.bet = 1 / (S.kB * S.Temp);
		str = 'undefined';    % initialize str
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'CHEB_DEGREE:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.npl = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'CHEFSI_OPTMZ:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.CheFSI_Optmz = C_param{1};
		fprintf('Neglecting option "%s", which is not supported in M-SPARC\n',str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'CHEFSI_BOUND_FLAG:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.chefsibound_flag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'FIX_RAND:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.FixRandSeed = C_param{1};
		fprintf('Neglecting option "%s", which is not supported in M-SPARC\n',str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'RHO_TRIGGER:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.rhoTrigger = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'NSTATES:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.Nev = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'NET_CHARGE:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.NetCharge = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		% error exit
		%msparc_error_exit(str);
	elseif (strcmp(str,'MAXIT_SCF:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MAXIT_SCF = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MINIT_SCF:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MINIT_SCF = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MAXIT_POISSON:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MAXIT_POISSON = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'RELAX_NITER:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.max_relax_it = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'RELAX_MAXDILAT:'))	
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);	
		S.max_dilatation = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line	
	elseif (strcmp(str,'TOL_RELAX_CELL:'))	
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);	
		S.TOL_RELAX_CELL = C_param{1};	
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line    
	elseif (strcmp(str,'ACCURACY:')) 
		C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
		temp = char(C_param{:});
		if (strcmp(temp,'minimal'))
			S.accuracy_level = 0;
		elseif (strcmp(temp,'low'))
			S.accuracy_level = 1;
		elseif (strcmp(temp,'medium'))
			S.accuracy_level = 2;
		elseif (strcmp(temp,'high'))
			S.accuracy_level = 3;
		elseif (strcmp(temp,'extreme'))
			S.accuracy_level = 4;
		else 
			error('\nCannot recognize accuracy level: "%s"\n',temp);
		end
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'SCF_FORCE_ACC:'))
		Flag_accuracy = Flag_accuracy + 1;
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.target_force_accuracy = C_param{1};
		%str = 'undefined';    % initialize str
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'SCF_ENERGY_ACC:')) 
		Flag_accuracy = Flag_accuracy + 1;
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.target_energy_accuracy = C_param{1};
		%str = 'undefined';    % initialize str
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'TOL_SCF:')) 
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.SCF_tol = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'TOL_POISSON:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.poisson_tol = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'TOL_RELAX:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.TOL_RELAX = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'TOL_LANCZOS:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.TOL_LANCZOS = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'TOL_PSEUDOCHARGE:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.pseudocharge_tol = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'TOL_KERKER:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		%S.kerker_tol = C_param{1};
		fprintf('WARNING: TOL_KERKER is obsolete, use TOL_PRECOND instead!\n')
		S.precond_tol = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'TOL_PRECOND:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.precond_tol = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PRECOND_KERKER_KTF:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.precond_kTF = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PRECOND_KERKER_THRESH:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.precond_thresh = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    elseif (strcmp(str,'PRECOND_KERKER_KTF_MAG:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.precond_kTF_mag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PRECOND_KERKER_THRESH_MAG:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.precond_thresh_mag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PRECOND_RESTA_Q0:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.precond_resta_q0 = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PRECOND_RESTA_RS:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.precond_resta_Rs = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PRECOND_FITPOW:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.precond_fitpow = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'REFERENCE_CUTOFF:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.rc_ref = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MIXING_VARIABLE:'))
		C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
		temp = char(C_param{:});
		if (strcmp(temp,'density'))
			S.MixingVariable = 0;
		elseif (strcmp(temp,'potential'))
			S.MixingVariable = 1;
		else
			fprintf('\nCannot recognize mixing variable: "%s"\n',temp);
			fprintf('Available options: "density", "potential" \n');
			error('exiting');
		end
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MIXING_PRECOND:'))
		C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
		temp = char(C_param{:});
		if (strcmp(temp,'none'))
			S.MixingPrecond = 0;
		elseif (strcmp(temp,'kerker'))
			S.MixingPrecond = 1;
		elseif (strcmp(temp,'resta'))
			S.MixingPrecond = 2;
		elseif (strcmp(temp,'truncated_kerker'))
			S.MixingPrecond = 3;
		else
			fprintf('\nCannot recognize mixing preconditioner: "%s"\n',temp);
			fprintf('Available options: "none", "kerker" \n');
			error('exiting');
		end
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    
    elseif (strcmp(str,'MIXING_PRECOND_MAG:'))
		C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
		temp = char(C_param{:});
		if (strcmp(temp,'none'))
			S.MixingPrecondMag = 0;
		elseif (strcmp(temp,'kerker'))
			S.MixingPrecondMag = 1;
		else
			fprintf('\nCannot recognize mixing preconditioner: "%s"\n',temp);
			fprintf('Available options: "none", "kerker" \n');
			error('exiting');
		end
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MIXING_HISTORY:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MixingHistory = C_param{1};
		%fscanf(fid1,'%f',S.MixingHistory);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MIXING_PARAMETER:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MixingParameter = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MIXING_PARAMETER_SIMPLE:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MixingParameterSimple = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    elseif (strcmp(str,'MIXING_PARAMETER_MAG:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MixingParameterMag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MIXING_PARAMETER_SIMPLE_MAG:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MixingParameterSimpleMag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PULAY_FREQUENCY:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.PulayFrequency = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PULAY_RESTART:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.PulayRestartFlag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'TWTIME:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.TWtime = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MD_FLAG:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MDFlag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MD_METHOD:'))
		C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MDMeth = char(C_param{:});
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MD_TIMESTEP:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MD_dt = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'MD_NSTEP:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MD_Nstep = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'ION_TEMP:'))
		Flag_ionT = Flag_ionT + 1;
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.ion_T = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'ION_TEMP_END:'))
		%Flag_ionT_end = Flag_ionT_end + 1;
		%C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		%S.thermos_Tf = C_param{1};
		msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'ION_ELEC_EQT:'))
		Flag_eqT = Flag_eqT + 1;
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.ion_elec_eqT = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'ION_VEL_DSTR:'))
		%C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		%S.ion_vel_dstr = C_param{1};
		msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'QMASS:'))
		%C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		%S.qmass = C_param{1};
		msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PRINT_MDOUT:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.PrintMDout = C_param{1};
		%fscanf(fid1,'%f',S.PrintMDout);
		%msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'RELAX_FLAG:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.RelaxFlag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'RELAX_METHOD:'))
		C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.RelaxMeth = char(C_param{:});
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'NLCG_SIGMA:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.NLCG_sigma = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'L_HISTORY:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.L_history = C_param{:};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'L_FINIT_STP:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.L_finit_stp = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'L_MAXMOV:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.L_maxmov = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'L_AUTOSCALE:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.L_autoscale = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'L_LINEOPT:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.L_lineopt = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'L_ICURV:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.L_icurv = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'FIRE_dt:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.FIRE_dt = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'FIRE_mass:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.FIRE_mass = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'FIRE_maxmov:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.FIRE_maxmov = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PRINT_RELAXOUT:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.PrintRelaxout = C_param{1};
		%msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'RESTART_FLAG:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.RestartFlag = C_param{1};
		%msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PRINT_RESTART:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.Printrestart = C_param{1};
		%msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PRINT_RESTART_FQ:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.Printrestart_fq = C_param{1};
		%msparc_neglect_warning(str);
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'EXCHANGE_CORRELATION:'))
		C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.XC = char(C_param{:});
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'CALC_STRESS:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.Calc_stress = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'CALC_PRES:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.Calc_pres = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'D3_FLAG:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.d3Flag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'D3_RTHR:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.d3Rthr = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'D3_CN_THR:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.d3Cn_thr = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'VDWDF_GEN_KERNEL:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.vdWDFKernelGenFlag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'NTYPES:'))
		fprintf('WARNING: NTYPES is no longer needed, skipping this option ...\n');
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'ATOMIC_MASS:'))
		fprintf('WARNING: ATOMIC_MASS is now moved to .ion file, skipping this option ...\n');
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'PSEUDOPOTENTIAL_LOCAL:'))
		fprintf('WARNING: PSEUDOPOTENTIAL_LOCAL is no longer needed, skipping this option ...\n');
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line% remove this, and this if cond
	elseif(strcmp(str,'PSEUDOPOTENTIAL_FILE:'))
		fprintf('WARNING: PSEUDOPOTENTIAL_FILE is now moved to .ion file, skipping this option ...\n');
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif(strcmp(str,'PRINT_FORCES:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.PrintForceFlag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif(strcmp(str,'PRINT_ATOMS:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.PrintAtomPosFlag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif(strcmp(str,'PRINT_VELS:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.PrintAtomVelFlag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line    
	elseif(strcmp(str,'PRINT_EIGEN:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.PrintEigenFlag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif(strcmp(str,'PRINT_DENSITY:'))
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.PrintElecDensFlag = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	elseif (strcmp(str,'OUTPUT_FILE:'))	  
		C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.filename_out = char(C_param{:});
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    elseif (strcmp(str,'MAXIT_FOCK:'))	  
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.MAXIT_FOCK = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    elseif (strcmp(str,'TOL_FOCK:'))	  
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.FOCK_TOL = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    elseif (strcmp(str,'TOL_SCF_INIT:'))	  
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.SCF_tol_init = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    elseif (strcmp(str,'EXX_METHOD:'))	  
		C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.ExxMethod = char(C_param{:});
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    elseif (strcmp(str,'ACE_FLAG:'))	  
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.ACEFlag = char(C_param{:});
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    elseif (strcmp(str,'EXX_ACE_VALENCE_STATES:'))	  
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.EXXACEVal_state = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    elseif (strcmp(str,'EXX_DOWNSAMPLING:'))
		C_param = textscan(fid1,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.exx_downsampling(1) = C_param{1};
		S.exx_downsampling(2) = C_param{2};
		S.exx_downsampling(3) = C_param{3};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    elseif (strcmp(str,'EXX_DIVERGENCE:'))
		C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.ExxDivMethod = char(C_param{:});
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    elseif (strcmp(str,'EXX_RANGE_FOCK:'))	  
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.hyb_range_fock = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
    elseif (strcmp(str,'EXX_RANGE_PBE:'))	  
		C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
		S.hyb_range_pbe = C_param{1};
		textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
	else 
		error('\nCannot recognize input variable identifier: "%s"\n',str);
		%fprintf('\nCannot recognize input flag in .inpt file: "%s"\n',str);
	end
end
	
fclose(fid1);

% check if the inpt options have conflict

if (S.Nx > 0 && S.Ny > 0 && S.Nz > 0) + (S.ecut > 0) + (S.mesh_spacing > 0) > 1
	error('"FD_GRID", "ECUT" and "MESH_SPACING" cannot be specified simultaneously!');
end

% check if non-orthogonal unit cell is used for Dirichlet BC
% if (S.BC == 1 && S.cell_typ == 2)
%     error('\nOnly orthogonal cells permitted for isolated clusters\n');
% end
n_Dirichlet = S.BCx + S.BCy + S.BCz;
if (S.BC >= 0) && (S.BCx >= 0 && S.BCy >= 0 && S.BCz >= 0)
	error('"BOUNDARY_CONDITION" and "BC" cannot be specified together!');
end

% % Boundary condition type: 1--isolated cluster; 2--periodic system
% if (S.cell_typ == 2) && (((S.BC == 1) || (S.BC == 4)) || (n_Dirichlet == 3 || n_Dirichlet == 2))
%     error('Non-orthogonal option allowed only for boundary conditions 2 & 3');
% end

% check if MD and Relaxation are turned on simultaneously
if (S.MDFlag ~= 0 && S.RelaxFlag ~= 0) 
	fprintf('\nStructural relaxations and MD cannot be turned on simultaneously!\n');
end

if(S.MDFlag == 1 && Flag_ionT == 0)
	error('\nIonic temperature must be specified for MD!\n');
end

if(S.MDFlag == 1 && S.ion_elec_eqT == 1)
	S.Temp = S.ion_T;
	S.bet = 1 / (S.kB * S.Temp);
end

% check CELL and LATVEC_SCALE
if Flag_cell == 1 && Flag_latvec_scale == 1
    error('\nCELL and LATVEC_SCALE cannot be specified simultaneously!\n');
end

% LACVEC_SCALE takes into account the length of the LATVEC's, so we'll scale the cell lengths
if Flag_latvec_scale == 1
    S.L1 = S.latvec_scale_x * norm(S.lat_vec(1,:));
    S.L2 = S.latvec_scale_y * norm(S.lat_vec(2,:));
    S.L3 = S.latvec_scale_z * norm(S.lat_vec(3,:));
end

if Flag_kptshift == 0
    S.kptshift = 0.5*(1-mod(S.nkpt,2));
end
end

function msparc_neglect_warning(str)
	fprintf('Neglecting option "%s", which is not supported in M-SPARC\n',str);
end

function msparc_error_exit(str)
	error('The option "%s" is not supported in M-SPARC\n',str);
end

