function S = md(S)
% @ brief   Function to perform molecular dynamics
% @ authors
%         Abhiraj Sharma <asharma424@gatech.edu>
%         Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%   Read the .restart file if MD is restarted
	if S.RestartFlag == 1
		S = RestartMD(S);
	end
	
	% set up MD output filenames
	S.mdfname = strcat(S.filename,'.aimd'); 
	if S.suffixNum > 0
		S.mdfname = sprintf('%s.aimd_%d',S.filename,S.suffixNum);
	end

	% Different MD methods
	if strcmp(S.MDMeth,'NVE')
		S = MD_NVE(S);
	else
		error('Currently only NVE is implemented in this code!');
	end
end




function S = MD_NVE(S)
% @ brief    Molecular dynamics of a canonical ensemble is performed
%            using velocity verlet algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% MD parameters
	MD_dt = S.MD_dt; % in femto-seconds
	MD_dt = MD_dt * S.fs2atu; % in atu
	kB = S.kB;

	% Atomic mass
	S.atomM = zeros(S.n_atm,1);
	count_typ = 1;
	count_typ_atms = 1;
	for JJ_a = 1:S.n_atm % loop over all the atoms
		S.atomM(JJ_a) = S.Atm(count_typ).Mass;
		if count_typ_atms == S.Atm(count_typ).n_atm_typ
			count_typ_atms = 1;
			count_typ = count_typ + 1;
		else
			count_typ_atms = count_typ_atms + 1;
		end
	end

	S.atomM = S.atomM * S.amu2au; % in au

	% Degree of freedom in the system
	dof = 3*S.n_atm - 3;

	% MD initialization
	if ~S.RestartFlag
		%x = transpose(S.lat_uvec) * S.Atoms';
		%x = reshape(x,[],1);
		x = coordinateTransformation(S, S.Atoms, 'noncart2cart_dis');
		x = reshape(x',[],1);
		% 1 - S.MD_velocity initialization (Maxwell-Boltzmann distr.)
		S.MD_vel = zeros(S.n_atm,3);
		rng('default');
		S.MD_vel(:,1) = sqrt(kB * S.ion_T./S.atomM) .* cos(2 * pi * rand(S.n_atm,1)) .* sqrt(-2.0 * log(rand(S.n_atm,1)));
		S.MD_vel(:,2) = sqrt(kB * S.ion_T./S.atomM) .* cos(2 * pi * rand(S.n_atm,1)) .* sqrt(-2.0 * log(rand(S.n_atm,1)));
		S.MD_vel(:,3) = sqrt(kB * S.ion_T./S.atomM) .* cos(2 * pi * rand(S.n_atm,1)) .* sqrt(-2.0 * log(rand(S.n_atm,1)));
		mvsum = sum(S.atomM' * S.MD_vel,1);
		
		% Remove translation
		S.MD_vel = S.MD_vel - repmat(mvsum/sum(S.atomM),S.n_atm,1);
		KE = 0.5 * sum(S.atomM' * (S.MD_vel .* S.MD_vel));

		% Rescale and reshape S.MD_velocities
		S.MD_vel = S.MD_vel * sqrt(dof * kB * S.ion_T/(2*KE));
	else
		x = reshape(S.Atoms',[],1);
	end

	% 2 - acceleration initialization
	[x,f,S] = electronicGroundStateAtomicForce(x,S); f = -f;
	S.ForceCount = S.ForceCount + 1;

	f = transpose(reshape(f,3,[]));
	acc = f./repmat(S.atomM,1,3);
	acc = reshape(acc',[],1);

	% Initial MD quantities
	KE = 0.5 * sum(S.atomM' * (S.MD_vel .* S.MD_vel));
	S.MD_vel = reshape(S.MD_vel',[],1);

	TE = (S.Etotal + KE)./S.n_atm;

	fprintf('Totalenergy/atom and temperature in MD step %d are: %.15f %.15f\n\n',1,TE, S.ion_T);

	 % Print the output in .aimd file
	if (S.PrintMDout == 1 && ~S.RestartFlag)
		PrintMD(S,transpose(reshape(x,3,[])),transpose(reshape(S.MD_vel,3,[])),TE);
	end
	S.Relax_iter = S.Relax_iter + 1;


	for i = 1:S.MD_Nstep+S.RestartFlag-1
		% Leapfrog step - I
		S.MD_vel = S.MD_vel + 0.5 * MD_dt * acc;
		x = x + MD_dt * S.MD_vel;

		% Calculate force
		[x,f,S] = electronicGroundStateAtomicForce(x,S); f = -f;
		S.ForceCount = S.ForceCount + 1;

		% Leapfrog step - II
		f = transpose(reshape(f,3,[]));
		acc = f./repmat(S.atomM,1,3);
		acc = reshape(acc',[],1);
		S.MD_vel = S.MD_vel + 0.5 * MD_dt * acc;
		S.MD_vel = transpose(reshape(S.MD_vel,3,[]));
		KE = 0.5 * sum(S.atomM' * (S.MD_vel .* S.MD_vel));
		S.MD_vel = reshape(S.MD_vel',[],1);
		TE = (S.Etotal + KE)/S.n_atm;
		S.ion_T = 2 * KE /(kB * dof);
		if(S.ion_elec_eqT == 1)
			S.Temp = S.ion_T;
			S.bet = 1 / (S.kB * S.Temp);  
		end
		fprintf('Totalenergy/atom and temperature in MD step %d are: %.15f %.15f\n\n',i+1,TE, S.ion_T);

		% print MD output in .aimd file
		if S.PrintMDout == 1
			PrintMD(S,transpose(reshape(x,3,[])),transpose(reshape(S.MD_vel,3,[])),TE);
		end
		% print relevant quantities for MD restart in .restart file 
		if (S.Printrestart == 1 && rem(S.Relax_iter, S.Printrestart_fq) == 0) 
			PrintRestartMD(S,transpose(reshape(x,3,[])),transpose(reshape(S.MD_vel,3,[])));
		end

		S.Relax_iter = S.Relax_iter + 1;
	end

	if (S.Printrestart == 1)
		S.Relax_iter = S.Relax_iter - 1;
		PrintRestartMD(S,transpose(reshape(x,3,[])),transpose(reshape(S.MD_vel,3,[])));
	end
end




function [] = PrintMD(S,atom_pos,atom_vel,TE)
% @ brief     Function to write MD output in .aimd file
%========================================================    
	if S.Relax_iter == 1
		fileID = fopen(S.mdfname,'w');
	else
		fileID = fopen(S.mdfname,'a');
	end
	fprintf(fileID,':MDSTEP: %d\n',S.Relax_iter);
	if (S.PrintAtomPosFlag == 1)
		fprintf(fileID,':R:\n');
		fprintf(fileID, '%18.10E %18.10E %18.10E\n', atom_pos');
	end
	if (S.PrintAtomVelFlag == 1)
		fprintf(fileID,':V:\n');
		fprintf(fileID, '%18.10E %18.10E %18.10E\n', atom_vel');
	end
	if (S.PrintForceFlag == 1)
		fprintf(fileID,':F:\n');
		fprintf(fileID, '%18.10E %18.10E %18.10E\n', S.force');
	end
	fprintf(fileID,':TEL: %18.10E\n',S.Temp);
	fprintf(fileID,':TIO: %18.10E\n',S.ion_T);
	fprintf(fileID,':TE: %18.10E\n',TE);
	fclose(fileID);
end




function [] = PrintRestartMD(S,atom_pos,atom_vel)
% @ brief     Function to write relevant quantities for MD restart in .restart file
%===================================================================================    
	fileID = fopen(S.restartfname,'w');
	fprintf(fileID,':MDSTEP: %d\n',S.Relax_iter);
	fprintf(fileID,':R:\n');
	fprintf(fileID, '%18.10E %18.10E %18.10E\n', atom_pos');
	fprintf(fileID,':V:\n');
	fprintf(fileID, '%18.10E %18.10E %18.10E\n', atom_vel');
	fprintf(fileID,':TEL: %18.10E\n',S.Temp);
	fprintf(fileID,':TIO: %18.10E\n',S.ion_T);
	fclose(fileID);
end



function S = RestartMD(S)
% @ brief        Function to read .restart file of MD
% ====================================================    
	fprintf(' Reading .restart file ...\n');
	fid = fopen(S.restartfname,'r');
	if (fid == -1) 
		error('\n Cannot open file "%s"\n',S.restartfname);
	end

	while(~feof(fid))
		C_inpt = textscan(fid,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
		str = char(C_inpt{:});
		if (strcmp(str,':MDSTEP:'))
			C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
			S.Relax_iter = C_param{1};
			textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		elseif (strcmp(str,':R:'))
			textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
			for i = 1:S.n_atm
				C_param = textscan(fid,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
				S.Atoms(i,:) = cell2mat(C_param);
			end
		elseif (strcmp(str,':V:'))
			S.MD_vel = zeros(S.n_atm,3);
			textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
			for i = 1:S.n_atm
				C_param = textscan(fid,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
				S.MD_vel(i,:) = cell2mat(C_param);
			end
		elseif (strcmp(str,':TEL:'))
			C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
			S.Temp = C_param{1};
			textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		elseif (strcmp(str,':TIO:'))
			C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
			S.ion_T = C_param{1};
			textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		end
	end
	S.bet = 1 / (S.kB * S.Temp);
	fclose(fid);
end
