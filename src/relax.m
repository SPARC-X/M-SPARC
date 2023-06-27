function S = relax(S)
% @ brief   Function to perform structural relaxation
% @ authors
%         Abhiraj Sharma <asharma424@gatech.edu>
%         Qimen Xu <qimenxu@gatech.edu>
%         Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%======================================================================

	% Read .restart file if relaxation is restarted
	if S.RestartFlag == 1
	  S = RestartRelax(S);
	end
	
	% set up relaxation output filename
	if S.RelaxFlag == 1
		S.relaxfname = strcat(S.filename,'.geopt'); 
		if S.suffixNum > 0
			S.relaxfname = sprintf('%s.geopt_%d',S.filename,S.suffixNum);
		end
	elseif S.RelaxFlag == 2
		S.relaxfname = strcat(S.filename,'.cellopt'); 
		if S.suffixNum > 0
			S.relaxfname = sprintf('%s.cellopt_%d',S.filename,S.suffixNum);
		end
	end
	
	% Different relaxation methods
	if S.RelaxFlag == 1
		if strcmp(S.RelaxMeth,'NLCG')
			S = NLCG(S);
		elseif strcmp(S.RelaxMeth,'LBFGS')
			S = LBFGS(S);
		elseif strcmp(S.RelaxMeth,'FIRE')
			S = FIRE(S);
		else
			error('Invalid relaxation method option! Current version only has NLCG, LBFGS, and FIRE algorithms.');
		end
	elseif S.RelaxFlag == 2
		S = cell_relax(S);
	end
end




function S = NLCG(S)
% @ brief     Function to perform structural relaxation using NON-LINEAR CONJUGATE GRADIENT method
% @ reference: "An introduction to the conjugate gradient method without the agonizing pain"
%===================================================================================================
	% Obtain cartesian coordinates and reshape into a vector
	if ~S.RestartFlag
		%x = transpose(S.lat_uvec) * S.Atoms';
		x = coordinateTransformation(S, S.Atoms, 'noncart2cart_dis');
		x = reshape(x',[],1);
	else
		x = reshape(S.Atoms',[],1);
	end
	
	imax = S.max_relax_it;
	jmax = 6;
	n = 30;
	tol1 = S.TOL_RELAX;
	tol2 = tol1*1e-2;
	sigma_0 = S.NLCG_sigma;
	i = 1;
	
	% Residual is same as atomic force (but the function returns dE/dR)
	[x,r,S] = electronicGroundStateAtomicForce(x,S); r = -r;
	S.ForceCount = S.ForceCount + 1;
	
	% Print the output in .geopt file
	if (S.PrintRelaxout == 1 && ~S.RestartFlag)
		PrintRelax(S,transpose(reshape(x,3,[])));
	end

	S.Relax_iter = S.Relax_iter + 1;
	k = S.Relax_iter;

	M = ones(size(x)); % TODO: Need to set REAL PRECONDITIONER
	s = r./M;
	delta_new = dot(r,s);
	% Search direction
	if ~S.RestartFlag
		S.NLCG_d = s;
	end
	% Error is defined as the supremum norm of the residual vector
	err = max(abs(r));

	while (i<imax && err>tol1)
		delta_d = max(abs(S.NLCG_d));
		[~,r_temp,S] = electronicGroundStateAtomicForce(x+sigma_0*S.NLCG_d,S);
		S.ForceCount = S.ForceCount + 1;
		eta_prev = dot(r_temp,S.NLCG_d);
		j = 0 ;
		alph = -sigma_0;
		r_temp = -r;
		% Find the step length (alph) using secant method
		while ((j < jmax) && (abs(alph)*delta_d > tol2))
			eta = dot(r_temp,S.NLCG_d);
			alph = alph * (eta / (eta_prev - eta));
			x = x + alph * S.NLCG_d;
			[x,r_temp,S] = electronicGroundStateAtomicForce(x,S);
			S.ForceCount = S.ForceCount + 1;
			eta_prev = eta;
			j = j + 1;
			fprintf('\n Secant step no. %d completed. \n',j);
		end
		delta_old = delta_new;
		r = -r_temp;        
		delta_mid = dot(r,s);
		M = ones(size(x)); %TODO: Need to set REAL PRECONDITIONER
		s = r./M;
		delta_new = dot(r,s);
		err = max(abs(r));
		bet = (delta_new-delta_mid)/delta_old;
		if (k==n || bet <=0)
			S.NLCG_d = s;
			k = 0;
		else
			S.NLCG_d = s + bet*S.NLCG_d;
		end

		% print the result in .geopt file
		if S.PrintRelaxout == 1
			PrintRelax(S,transpose(reshape(x,3,[])));
		end
		%print relevant quanties for restart in .restart file
		if (S.Printrestart == 1 && rem(S.Relax_iter, S.Printrestart_fq) == 0) 
			PrintRestart(S,transpose(reshape(x,3,[])));
		end

		k = k + 1;
		i = i + 1;
		S.Relax_iter = S.Relax_iter + 1;
		fprintf('\n Outer CG iteration no. %d completed. \n',i);
		fprintf('\n error in CG relaxation = %.9f\n',delta_new);
	end

	if (S.Printrestart == 1)
		S.Relax_iter = S.Relax_iter - 1;
		PrintRestart(S,transpose(reshape(x,3,[])));
	end
	fprintf('\n Total no. of outer CG iterations = %d. ',i);

end




function S = LBFGS(S) 
% @ brief  Function to perform relaxation using LIMITED-MEMORY BFGS 
%          (Broyden–Fletcher–Goldfarb–Shanno) method
%@ reference: "Based on the implementation in VTST(VASP Transition State Tools)"
%================================================================================
	imax = S.max_relax_it;
	lbfgs_tol = S.TOL_RELAX;
	N = 3 * S.n_atm;

	% User input
	M = S.L_history;
	FINIT_STP = S.L_finit_stp;
	MAXMOV = S.L_maxmov;
	AUTOSCALE = S.L_autoscale;
	LINEOPT = S.L_lineopt; % needed only if AUTOSCALE = 0
	ICURV = S.L_icurv; % needed only if AUTOSCALE = 0;

	DAMPING = 2.0; % Don't change
	if ~S.RestartFlag
		S.LBFGS_ISFD = 1; % Don't change
		S.LBFGS_ISRESET = 1; % Don't change
		S.LBFGS_STEP = 0; % Don't change

		S.LBFGS_DELTAX = zeros(N,M);
		S.LBFGS_DELTAG = zeros(N,M);
		S.LBFGS_IYS = zeros(M,1);
		S.LBFGS_FOLD = zeros(N,1);
		S.LBFGS_ATOMDISP = zeros(N,1);
	end
	ALPHA = zeros(M,1);
	
	% Obtain atomic positions in cartesian coordinates in a vector form
	if ~S.RestartFlag
		%X = transpose(S.lat_uvec) * S.Atoms';
		X = coordinateTransformation(S, S.Atoms, 'noncart2cart_dis');
		X = reshape(X',[],1);
	else
		X = reshape(S.Atoms',[],1);
	end 

	% Residual (F) is equal to atomic force (function returns dE/dR)
	[X,G,S] = electronicGroundStateAtomicForce(X,S); F = -G;
	S.ForceCount = S.ForceCount + 1;

	% Print the relaxation output
	if (S.PrintRelaxout == 1 && ~S.RestartFlag)
		PrintRelax(S,transpose(reshape(X,3,[])));
	end

	S.Relax_iter = S.Relax_iter + 1;   
	err = max(abs(F));
	i = 1;

	while(i < imax && err > lbfgs_tol)
		MAXMOV_FLAG = 0;
		if(AUTOSCALE == 1)
			if(S.LBFGS_ISFD == 1)
				DIR = F/norm(F);
				XOLD = X;
				S.LBFGS_FOLD = F;
				X = X + DIR * FINIT_STP;
				S.LBFGS_ISFD = 0;
			else
				% Compute inverse curvature 
				DIR = S.LBFGS_ATOMDISP;
				%disp(DIR);
				DIR_NRM = norm(DIR);
				%disp(DIR_NRM);
				FP1 = dot(S.LBFGS_FOLD, DIR)/DIR_NRM;
				FP2 = dot(F, DIR)/DIR_NRM;
				CURV = (FP1 - FP2)/DIR_NRM;
				%disp(CURV);
				ICURV = 1/CURV/DAMPING; % must be positive-definite
				if(ICURV < 0)
					S.LBFGS_ISRESET = 1;
					MAXMOV_FLAG = 1;
				end
				if(S.LBFGS_ISRESET == 1)
					S.LBFGS_STEP = 0;
					S.LBFGS_ISRESET = 0;
				else
					S.LBFGS_STEP = S.LBFGS_STEP + 1;
					if(S.LBFGS_STEP <= M)
						S.LBFGS_DELTAX(:,S.LBFGS_STEP) = S.LBFGS_ATOMDISP;
						S.LBFGS_DELTAG(:,S.LBFGS_STEP) = -(F - S.LBFGS_FOLD);
						S.LBFGS_IYS(S.LBFGS_STEP) = 1/dot(S.LBFGS_DELTAX(:,S.LBFGS_STEP),S.LBFGS_DELTAG(:,S.LBFGS_STEP));
					else
						S.LBFGS_DELTAX(:,1:M-1) = S.LBFGS_DELTAX(:,2:M);
						S.LBFGS_DELTAG(:,1:M-1) = S.LBFGS_DELTAG(:,2:M);
						S.LBFGS_IYS(1:M-1) = S.LBFGS_IYS(2:M);
						S.LBFGS_DELTAX(:,M) = S.LBFGS_ATOMDISP;
						S.LBFGS_DELTAG(:,M) = -(F - S.LBFGS_FOLD);
						S.LBFGS_IYS(M) = 1/dot(S.LBFGS_DELTAX(:,M),S.LBFGS_DELTAG(:,M));
					end
				end
				XOLD = X;
				S.LBFGS_FOLD = F;
				% Compute H0*g
				if(S.LBFGS_STEP <= M)
					BOUND = S.LBFGS_STEP;
				else
					BOUND = M;
				end
				DIR = -F;

				for IM = 1:BOUND
					JM =  BOUND + 1 - IM;
					ALPHA(JM) = S.LBFGS_IYS(JM) * dot(S.LBFGS_DELTAX(:,JM),DIR);
					DIR =  DIR - ALPHA(JM) * S.LBFGS_DELTAG(:,JM);
				end
				DIR = ICURV * DIR;
				for IM = 1:BOUND
					BETA = S.LBFGS_IYS(IM) * dot(S.LBFGS_DELTAG(:,IM), DIR);
					DIR = DIR + S.LBFGS_DELTAX(:,IM) * (ALPHA(IM) - BETA);
				end
				DIR = -DIR;
				%disp(DIR);
				STP_SIZ = norm(DIR);
				%disp(STP_SIZ);
				if(STP_SIZ > MAXMOV)
					S.LBFGS_ISRESET = 1;
					STP_SIZ = MAXMOV;
					DIR = STP_SIZ * F/norm(F);
				end
				if(MAXMOV_FLAG == 1)
					DIR = F/norm(F);
					X = X + MAXMOV * DIR;
					MAXMOV_FLAG = 0;
				else
					X = X + DIR;
				end
			   % disp(X);
			end
		else
			if(S.LBFGS_ISFD == 1)
				S.LBFGS_ISFD = 0;
				OPTFLAG = 1;
				a1 = abs(dot(F,S.LBFGS_FOLD));
				a2 = dot(S.LBFGS_FOLD,S.LBFGS_FOLD);
				if(a1 > 0.5 * a2 || a2 == 0)
					S.LBFGS_ISRESET = 1;
				end
				if(LINEOPT == 0)
					S.LBFGS_ISRESET = 0;
				end
				if(a2 == 0)
					S.LBFGS_ISRESET = 1;
				end
				if(S.LBFGS_ISRESET == 1)
					S.LBFGS_STEP = 0;
					S.LBFGS_ISRESET = 0;
				else
					S.LBFGS_STEP = S.LBFGS_STEP + 1;
					if(S.LBFGS_STEP <= M)
						S.LBFGS_DELTAX(:,S.LBFGS_STEP) = S.LBFGS_ATOMDISP;
						S.LBFGS_DELTAG(:,S.LBFGS_STEP) = -(F - S.LBFGS_FOLD);
						S.LBFGS_IYS(S.LBFGS_STEP) = 1/dot(S.LBFGS_DELTAX(:,S.LBFGS_STEP),S.LBFGS_DELTAG(:,S.LBFGS_STEP));
					else
						S.LBFGS_DELTAX(:,1:M-1) = S.LBFGS_DELTAX(:,2:M);
						S.LBFGS_DELTAG(:,1:M-1) = S.LBFGS_DELTAG(:,2:M);
						S.LBFGS_IYS(1:M-1) = S.LBFGS_IYS(2:M);
						S.LBFGS_DELTAX(:,M) = S.LBFGS_ATOMDISP;
						S.LBFGS_DELTAG(:,M) = -(F - S.LBFGS_FOLD);
						S.LBFGS_IYS(M) = 1/dot(S.LBFGS_DELTAX(:,M),S.LBFGS_DELTAG(:,M));
					end
				end
				XOLD = X;
				S.LBFGS_FOLD = F;
				% Compute H0*g
				if(S.LBFGS_STEP <= M)
					BOUND = S.LBFGS_STEP;
				else
					BOUND = M;
				end
				DIR = -F;
				for IM = 1:BOUND
					JM =  BOUND + 1 - IM;
					ALPHA(JM) = S.LBFGS_IYS(JM) * dot(S.LBFGS_DELTAX(:,JM),DIR);
					DIR =  DIR - ALPHA(JM) * S.LBFGS_DELTAG(:,JM);
				end
				DIR = ICURV * DIR;
				for IM = 1:BOUND
					BETA = S.LBFGS_IYS(IM) * dot(S.LBFGS_DELTAG(:,IM), DIR);
					DIR = DIR + S.LBFGS_DELTAX(:,IM) * (ALPHA(IM) - BETA);
				end
				DIR = -DIR;
				if(LINEOPT == 1)
					DIR = DIR/norm(DIR);
					X = X + DIR * FINIT_STP;
				else
					STP_SIZ = norm(DIR);
					if(STP_SIZ > MAXMOV)
						STP_SIZ = MAXMOV;
						DIR = STP_SIZ * DIR /norm(DIR);
					end
					X = X + DIR;
					S.LBFGS_ISFD = 1;
					OPTFLAG = 0;
				end
			else
				S.LBFGS_ISFD = 1;
				OPTFLAG = 0;
				FP1 = dot(S.LBFGS_FOLD,DIR);
				FP2 = dot(F,DIR);
				CURV = (FP1 - FP2)/FINIT_STP;
				if(CURV < 0)
					STP_SIZ = MAXMOV;
				else
					FAVG = 0.5 * (FP1 + FP2);
					STP_SIZ = FAVG/CURV;
					if(abs(STP_SIZ) > MAXMOV)
						if(STP_SIZ >= 0)
							STP_SIZ = abs(MAXMOV) - abs(FINIT_STP);
						else
							STP_SIZ = -abs(MAXMOV) + abs(FINIT_STP);
						end
					else
						STP_SIZ =  STP_SIZ - 0.5 * FINIT_STP;
					end
				end
				X = X + DIR * STP_SIZ;
			end
		end
		%X(1:3) = Y(1:3);
		S.LBFGS_ATOMDISP = X - XOLD;
		[X,G,S] = electronicGroundStateAtomicForce(X,S); F = -G;
		S.ForceCount = S.ForceCount + 1;
		err = max(abs(F));

		% print the relaxed configuration in .geopt file
		if S.PrintRelaxout == 1
			PrintRelax(S,transpose(reshape(X,3,[])));
		end
		% print the relevant quantities needed for restart in .restart file
		if (S.Printrestart == 1 && rem(S.Relax_iter, S.Printrestart_fq) == 0) 
			PrintRestart(S,transpose(reshape(X,3,[])));
		end

		i = i + 1;
		S.Relax_iter = S.Relax_iter + 1;
	end
	
	if (S.Printrestart == 1)
		S.Relax_iter = S.Relax_iter - 1;
		PrintRestart(S,transpose(reshape(X,3,[])));
	end
end





function S = FIRE(S)
% @ brief      Function to perform relaxation using  FAST INERTIAL RELAXATION ENGINE method
% @ reference:  "Structural Relaxation Made Simple (Bitzek et. al., 2006)"
%==========================================================================================
	% Obtain atomic positions in cartesian coordinates and as a vector
	if ~S.RestartFlag
		%x_var = transpose(S.lat_uvec) * S.Atoms';
		x_var = coordinateTransformation(S, S.Atoms, 'noncart2cart_dis');
		x_var = reshape(x_var',[],1);
	else
		x_var = reshape(S.Atoms',[],1);
	end

	% FIRE related variables : user modifiable
	S.FIRE_dt = S.FIRE_dt * S.fs2atu;
	FIRE_Mass = S.FIRE_mass * S.amu2au;
	FIRE_Tol = S.TOL_RELAX;
	Max_FIRE_Iter = S.max_relax_it;
	dx_Max = S.FIRE_maxmov;

	% FIRE related variables : internal
	fDec = 0.5;
	fInc = 1.1;
	nMin = 5;
	alphaStart = 0.1;
	fAlpha = 0.99;
	dt_Max = 10.0 * S.FIRE_dt;

	% Do a static calculation first
	%fprintf('\n\n Starting single shot static calculation: \n\n')
	[x_var,force_fun,S] = electronicGroundStateAtomicForce(x_var,S); Force_var = -force_fun;
	S.ForceCount = S.ForceCount + 1;

	% print atomic configuration in .geopt file
	if (S.PrintRelaxout == 1 && ~S.RestartFlag)
		PrintRelax(S,transpose(reshape(x_var,3,[])));
	end

	S.Relax_iter = S.Relax_iter + 1;
	
	if ~S.RestartFlag
		S.FIRE_resetIter = 0;
		S.FIRE_alpha = alphaStart;
		S.FIRE_dtNow = S.FIRE_dt;
		S.FIRE_vel = 0.0 * Force_var; % Set zero initial velocities
	end    

	err = max(abs(Force_var));
	FIRE_iter = S.Relax_iter;
	iter = 1;
	
	while (iter < Max_FIRE_Iter && err > FIRE_Tol)

		% Verlet step
		S.FIRE_vel = S.FIRE_vel + 0.5 * S.FIRE_dtNow * (Force_var/FIRE_Mass);
		delta_x = S.FIRE_vel * S.FIRE_dtNow;

		% Apply a "trust region"    
		ind = (delta_x > dx_Max);
		delta_x(ind) = dx_Max;

		ind = (delta_x < -dx_Max);
		delta_x(ind) = -dx_Max;

		x_var =  x_var + delta_x;
		[x_var,force_fun,S] = electronicGroundStateAtomicForce(x_var,S); Force_var = -force_fun;
		S.ForceCount = S.ForceCount + 1;

		S.FIRE_vel = S.FIRE_vel + 0.5 * S.FIRE_dtNow * (Force_var/FIRE_Mass);

		% Check power
		P = dot(Force_var, S.FIRE_vel);

		if(P < 0.0)

			S.FIRE_vel = 0.0 * S.FIRE_vel;  % Reset velosity to 0
			S.FIRE_resetIter = FIRE_iter;          % S.FIRE_resetIter <-- iter #
			S.FIRE_dtNow = S.FIRE_dtNow*fDec;   % decrease dt
			S.FIRE_alpha = alphaStart;         % reset alpha to alpha_start

		elseif ((P >= 0.0) && (( FIRE_iter - S.FIRE_resetIter) > nMin))

			S.FIRE_dtNow = min(S.FIRE_dtNow*fInc, dt_Max);   % Update dt
			S.FIRE_alpha = fAlpha*S.FIRE_alpha;        % update alpha

		end

		% Adjust velocities
		hatF = Force_var / norm(Force_var);   % the unit vector
		S.FIRE_vel = (1.0 - S.FIRE_alpha)*S.FIRE_vel + (S.FIRE_alpha * norm(S.FIRE_vel)) * hatF ;

		% print atomic configuration in .geopt file
		if S.PrintRelaxout == 1
			PrintRelax(S,transpose(reshape(x_var,3,[])));
		end
		% print relevant quantities needed for restart in .restart file
		if (S.Printrestart == 1 && rem(S.Relax_iter, S.Printrestart_fq) == 0) 
			PrintRestart(S,transpose(reshape(x_var,3,[])));
		end   

		err = max(abs(Force_var)); 

		FIRE_iter = FIRE_iter + 1;
		iter = iter + 1;  
		S.Relax_iter = S.Relax_iter + 1;     
	end

	if (S.Printrestart == 1)
		S.Relax_iter = S.Relax_iter - 1;
		PrintRestart(S,transpose(reshape(x_var,3,[])));
	end

end




function [] = PrintRelax(S,atom_pos)
% @ brief   Function to write relaxation output in .geopt file
%================================================================     
	if S.Relax_iter == 1
		fileID = fopen(S.relaxfname,'w');
	else
		fileID = fopen(S.relaxfname,'a');
	end
	fprintf(fileID,':RELAXSTEP: %d\n',S.Relax_iter);
	if (S.PrintAtomPosFlag == 1)
		fprintf(fileID,':R:\n');
		fprintf(fileID, '%18.10E %18.10E %18.10E\n', atom_pos');
	end
	if (S.PrintForceFlag == 1)
		fprintf(fileID,':F:\n');
		fprintf(fileID, '%18.10E %18.10E %18.10E\n', S.force');
	end
	fclose(fileID);
end




function [] = PrintRestart(S,atom_pos)
% brief     Function to write relevant quantities needed for restart of relaxation
%==================================================================================    
	fileID = fopen(S.restartfname,'w');
	fprintf(fileID,':RELAXSTEP: %d\n',S.Relax_iter);
	fprintf(fileID,':R:\n');
	fprintf(fileID, '%18.10E %18.10E %18.10E\n', atom_pos');
	
	if (strcmp(S.RelaxMeth,'NLCG'))
		% Print search direction
		fprintf(fileID,':D:\n');
		fprintf(fileID,'%18.10E\n', S.NLCG_d);   
	elseif (strcmp(S.RelaxMeth,'LBFGS'))
		fprintf(fileID,':ISFD: %d\n', S.LBFGS_ISFD);
		fprintf(fileID,':ISRESET: %d\n', S.LBFGS_ISRESET);
		fprintf(fileID,':STEP: %d\n', S.LBFGS_STEP);
		fprintf(fileID,':DX:\n');
		DX = reshape(S.LBFGS_DELTAX,[],1);
		fprintf(fileID,'%18.10E\n', DX);
		fprintf(fileID,':DG:\n');
		DG = reshape(S.LBFGS_DELTAG,[],1);
		fprintf(fileID,'%18.10E\n', DG);
		fprintf(fileID,':IYS:\n');
		fprintf(fileID,'%18.10E\n', S.LBFGS_IYS);
		fprintf(fileID,':FOLD:\n');
		fprintf(fileID,'%18.10E\n', S.LBFGS_FOLD);
		fprintf(fileID,':RDISP:\n');
		fprintf(fileID,'%18.10E\n', S.LBFGS_ATOMDISP);
	elseif(strcmp(S.RelaxMeth,'FIRE'))
		fprintf(fileID,':FIRE_alpha: %18.10E\n', S.FIRE_alpha);
		fprintf(fileID,':FIRE_dtNow: %18.10E\n', S.FIRE_dtNow);
		fprintf(fileID,':FIRE_resetIter: %d\n', S.FIRE_resetIter); 
		fprintf(fileID,':FIRE_V:\n');
		fprintf(fileID,'%18.10E\n', S.FIRE_vel);
	end

	fclose(fileID);
end




function S = RestartRelax(S)
% @ brief      Function to read the .restart file of relaxation
%================================================================    
	fprintf(' Reading .restart file ...\n');
	fid = fopen(S.restartfname,'r');
	if (fid == -1) 
		error('\n Cannot open file "%s"\n',S.restartfname);
	end

	while(~feof(fid))
		C_inpt = textscan(fid,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
		str = char(C_inpt{:});
		if (strcmp(str,':RELAXSTEP:'))
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
		elseif (strcmp(str,':D:'))
			S.NLCG_d = zeros(3*S.n_atm,1);
			textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
			for i = 1:3*S.n_atm
				C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
				S.NLCG_d(i) = C_param{1};
			end
		elseif (strcmp(str,':ISFD:'))
				C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
				S.LBFGS_ISFD = C_param{1};
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		elseif (strcmp(str,':ISRESET:'))
				C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
				S.LBFGS_ISRESET = C_param{1};
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		elseif (strcmp(str,':STEP:'))
				C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
				S.LBFGS_STEP = C_param{1};
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line                          
		elseif (strcmp(str,':DX:'))
				S.LBFGS_DELTAX = zeros(3*S.n_atm,S.L_history);
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
				for j = 1:S.L_history
					for i = 1:3*S.n_atm
						C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
						textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
						S.LBFGS_DELTAX(i,j) = C_param{1};
					end
				end 
		elseif (strcmp(str,':DG:'))
				S.LBFGS_DELTAG = zeros(3*S.n_atm,S.L_history);
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
				for j = 1:S.L_history
					for i = 1:3*S.n_atm
						C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
						textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
						S.LBFGS_DELTAG(i,j) = C_param{1};
					end
				end
		elseif (strcmp(str,':IYS:'))
				S.LBFGS_IYS = zeros(S.L_history,1);
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
				for i = 1:S.L_history
					C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
					textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
					S.LBFGS_IYS(i) = C_param{1};
				end
		elseif (strcmp(str,':FOLD:'))
				S.LBFGS_FOLD = zeros(3*S.n_atm,1);
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
				for i = 1:3*S.n_atm
					C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
					textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
					S.LBFGS_FOLD(i) = C_param{1};
				end
		elseif (strcmp(str,':RDISP:'))
				S.LBFGS_ATOMDISP = zeros(3*S.n_atm,1);
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
				for i = 1:3*S.n_atm
					C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
					textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
					S.LBFGS_ATOMDISP(i) = C_param{1};
				end
		elseif (strcmp(str,':FIRE_alpha:'))
				C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
				S.FIRE_alpha = C_param{1};
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		elseif (strcmp(str,':FIRE_dtNow:'))
				C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
				S.FIRE_dtNow = C_param{1};
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		elseif (strcmp(str,':FIRE_resetIter:'))
				C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
				S.FIRE_resetIter = C_param{1};
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
		elseif (strcmp(str,':FIRE_V:'))
				S.FIRE_vel = zeros(3*S.n_atm,1);
				textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
				for i = 1:3*S.n_atm
					C_param = textscan(fid,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
					textscan(fid,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0); % skip current line
					S.FIRE_vel(i) = C_param{1};
				end                                
		end
	end

	fclose(fid);
end















