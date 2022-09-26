function [atompos_next,AtmForce,S] = electronicGroundStateAtomicForce(atom_pos,S)
% @brief    electronicGroundStateAtomicForce(atom_pos,S) calculates the electronic ground
%           state atomic force corresponding to the given atom positions.
%
% @param atom_pos       Given atom positions (stretched into a (3*n_atoms)-by-1 vector)
% @param atompos_next   Atom position for the next relaxation/MD step (in cartesian coordinates)
% @param AtmForce       Atomic force as a 3*N x 1 vector 
% @authors	Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%==============================================================================================

tic_relax = tic;
fprintf('\n');
fprintf(' ###############################################################\n');
fprintf(' Relaxation step number: %d \n', S.Relax_iter);

if (S.RelaxFlag || S.MDFlag)
	% Reshaping atom positions
	atom_pos = transpose(reshape(atom_pos,3,[]));
	% update un-mapped atomic positions, for charge extrapolation
	S.atom_pos_tp1 = S.atom_pos_tp1 + (atom_pos - S.Atoms_old); 
	% S.dV_tp1 = S.dV;
	% Convert atom coordinates from cartesian to cell coordinates
	S.Atoms = coordinateTransformation(S, atom_pos, 'cart2noncart_dis');
end

% Reference_cutoff checks
nndis = calculate_min_distance(S);
if S.rc_ref > 0.5*nndis
    fprintf("\n WARNING: REFERENCE _CUFOFF (%.6f Bohr) > 1/2 nn (nearest neighbor) distance (%.6f Bohr) in SCF#%d\n",...
            S.rc_ref, 0.5*nndis, S.Relax_iter);
end
if S.rc_ref < S.dx || S.rc_ref < S.dy || S.rc_ref < S.dz
    fprintf("\n WARNIG: REFERENCE _CUFOFF (%.6f Bohr) < MESH_SPACING (dx %.6f Bohr, dy %.6f Bohr, dz %.6f Bohr) in SCF#%d\n",...
            S.rc_ref, S.dx, S.dy, S.dz, S.Relax_iter);
end
    
% Check position of atom near the boundary and apply wraparound in case of PBC
S = check_atomlocation(S);

% perform charge extrapolation (before rho_at is updated)
if (S.RelaxFlag == 1 || S.MDFlag)
	S = chargeExtrapolation(S);
end

% Pseudocharge (and reference), sum atomic charge density, self energy 
% (and reference), electrostatic correction 
% S.b,S.b_ref,S.Eself,S.Eself_ref,S.rho_at,S.E_corr,S.V_c, S.NegCharge,
% S.PosCharge, S.NetCharge
t_calc_b = tic;

S = calculate_b_guessRho_Eself(S);

fprintf(' Time for b calculation: %.3f seconds.\n',toc(t_calc_b));

% set up guess electron density (guess rho)
S = initElectrondensity(S);

% Calculate nonlocal projectors	
S.Atom = calculate_nloc_projector(S);

% Self-consistent Field (SCF) method
S = scf(S);
	
%save('rho.mat','-struct','S','rho');
S.S_Debug.relax(S.Relax_iter).Eself = S.Eself;
S.S_Debug.relax(S.Relax_iter).Eself_ref = S.Eself_ref;
S.S_Debug.relax(S.Relax_iter).E_corr = S.E_corr;

if abs(1-S.occ(1))>1e-6 || abs(S.occ(end))>1e-6
	fprintf('[\b WARNING: No. of states is not enough!]\b \n');
	S.S_Debug.relax(S.Relax_iter).occ_check = 1; % 1 means not satisfied
else
	S.S_Debug.relax(S.Relax_iter).occ_check = 0;
end

if S.d3Flag == 1
	if ispc % windows
		addpath('vdW\d3\');
	else % max/linux
		addpath('vdW/d3/');
	end
	S = d3EnergyGradient(S);
    S.Etotal = S.Etotal + S.d3Energy;
end
	
% Etotal = evaluateTotalEnergy(EigVal,occ,rho,S.b,phi,Vxc,S.W,S.bet,S.Eself,S.E_corr,1) ;

fprintf('\n');
fprintf(' **********************************************************\n');
fprintf(' *          Energy per unit cell = %11.9f Ha.       *\n', S.Etotal);
fprintf(' *          Energy per atom = %11.9f Ha.            *\n', S.Etotal / S.n_atm);
fprintf(' **********************************************************\n');

% write to output file
outfname = S.outfname;
fileID = fopen(outfname,'a');
if (fileID == -1) 
	error('\n Cannot open file "%s"\n',outfname);
end 
fprintf(fileID,'====================================================================\n');
fprintf(fileID,'                                Energy                              \n');
fprintf(fileID,'====================================================================\n');
fprintf(fileID,'Free energy per atom               :%18.10E (Ha/atom)\n', S.Etotal / S.n_atm);
fprintf(fileID,'Total free energy                  :%18.10E (Ha)\n', S.Etotal);
fprintf(fileID,'Band structure energy              :%18.10E (Ha)\n', S.Eband);
fprintf(fileID,'Exchange correlation energy        :%18.10E (Ha)\n', S.Exc);
fprintf(fileID,'Self and correction energy         :%18.10E (Ha)\n', S.E_corr-S.Eself);
fprintf(fileID,'Entropy*kb*T                       :%18.10E (Ha)\n', S.Eent);
fprintf(fileID,'Fermi level                        :%18.10E (Ha)\n', S.lambda_f);
if S.d3Flag == 1
	fprintf(fileID,'DFT-D3 correction                  :%18.10E (Ha)\n', S.d3Energy);
end
if (S.vdWDFFlag == 1) || (S.vdWDFFlag == 2)
    fprintf(fileID,'vdWDF energy                       :%18.10E (Ha)\n', S.vdWenergy);
end
if S.nspin ~= 1
	fprintf(fileID,'Net Magnetization                  :%18.10E \n', S.netM);
end
fclose(fileID);


% Atomic force calculation
tic_force = tic;
S.force = atomicForce(S);
%S.force = zeros(S.n_atm,3);
force_mat = S.force;
sz_fmat = size(force_mat);
force_corr = sum(force_mat, 1) / sz_fmat(1);
if S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5
	force_corr(1) = 0.0; % x - component
	force_corr(2) = 0.0; % y - component
end
S.force = force_mat - repmat(force_corr, sz_fmat(1),1);

% Apply constraint on atoms by making forces on fixed atoms zero
if S.RelaxFlag || S.MDFlag
	S.force = S.force .* S.mvAtmConstraint;
end

fprintf(' ***********************************************************\n');
fprintf(' *                      Atomic Force                       *\n');
fprintf(' ***********************************************************\n');
fprintf(' Drift free forces (Ha/Bohr):\n');
S.S_Debug.relax(S.Relax_iter).force = S.force;
disp(S.force);
S.abs_force = sqrt(sum(abs(S.force).^2,2));
fprintf(' Max magnitude of forces (Ha/Bohr):');
disp(max(S.abs_force));

t_force = toc(tic_force);
fprintf('\n Time for calculating forces: %f s.\n', t_force);

% write forces into .static file if required
if (S.PrintForceFlag == 1 && S.MDFlag == 0 && S.RelaxFlag == 0) 
	staticfname = S.staticfname;
	fileID = fopen(staticfname,'a');
	if (fileID == -1) 
		error('\n Cannot open file "%s"\n',staticfname);
	end
	fprintf(fileID, 'Atomic forces (Ha/Bohr):\n');
	fprintf(fileID, '%18.10f %18.10f %18.10f\n', S.force');
	fclose(fileID);
end

% write force magnitude and timing info. to .out file
outfname = S.outfname;
fileID = fopen(outfname,'a');
if (fileID == -1) 
	error('\n Cannot open file "%s"\n',outfname);
end

force_magnitude = sqrt(sum(S.force .* S.force, 2));
avgF = sum(force_magnitude) / length(force_magnitude);
maxF = max(force_magnitude);

fprintf(fileID,'Average force                      :%18.10E (Ha/Bohr)\n',avgF);
fprintf(fileID,'Maximum force                      :%18.10E (Ha/Bohr)\n',maxF);
fprintf(fileID,'Time for force calculation         :  %.3f (sec)\n',t_force);
fclose(fileID);


% Assign output of function
AtmForce = reshape(S.force',[],1);
AtmForce = - AtmForce;


% Perform stress and pressure calculation
if S.Calc_stress
	t1 = tic;
	Stress = evaluateStress(S);
	sigma_inGPa = Stress*2.94210119*(10^4);
	pressure = -trace(sigma_inGPa)/3.0;
	S.Stress = sigma_inGPa; % store stress tensor	
	S.Pressure = pressure;  % store pressure
	fprintf('\n[\b"Stress in GPa"\n\n\n]\b');
	disp(sigma_inGPa);
	fprintf('\n Time for calculating stress: %f s.\n\n', toc(t1));

	outfname = S.outfname;
	fileID = fopen(outfname,'a');
	if (fileID == -1) 
		error('\n Cannot open file "%s"\n',outfname);
	end
	if (S.MDFlag == 0 && S.RelaxFlag == 0)
		fprintf(fileID, 'Stress (GPa):\n');
		fprintf(fileID,'%18.10f %18.10f %18.10f\n', sigma_inGPa');
	end
	fprintf(fileID,'Pressure                           :%18.10E (GPa)\n',pressure);
	fprintf(fileID,'Maximum stress                     :%18.10E (GPa)\n',max(max(abs(sigma_inGPa))));
	fprintf(fileID,'Time for stress calculation        :  %.3f (sec)\n',toc(t1));
	
	fclose(fileID);        
end

if (S.Calc_pres && ~S.Calc_stress)
	t1 = tic;
	Pressure = evaluatePressure(S);
	PinGPa = Pressure*2.94210119*(10^4);
	S.Pressure = PinGPa;  % store pressure
	fprintf('\n[\b"Pressure %f Ha/Bohr^3 %f GPa"\n\n\n]\b',Pressure,PinGPa);
	fprintf('\n Time for calculating pressure: %f s.\n\n', toc(t1));

	outfname = S.outfname;
	fileID = fopen(outfname,'a');
	if (fileID == -1) 
		error('\n Cannot open file "%s"\n',outfname);
	end
	fprintf(fileID,'Pressure                           :%18.10E (GPa)\n',PinGPa);
	fprintf(fileID,'Time for pressure calculation      :  %.3f (sec)\n',toc(t1));
		
	fclose(fileID);
end 

% Updating relaxation iteration number
S.S_Debug.relax(S.Relax_iter).relax_time = toc(tic_relax);
if (S.RelaxFlag || S.MDFlag)
	fprintf('\n Relaxation step number %d completed in %f s.\n',S.Relax_iter, S.S_Debug.relax(S.Relax_iter).relax_time);

	outfname = S.outfname;
	fileID = fopen(outfname,'a');
	if (fileID == -1) 
		error('\n Cannot open file "%s"\n',outfname);
	end
	fprintf(fileID,'Relax time                         :  %.3f (sec)\n',toc(tic_relax));
	fclose(fileID);

end
fprintf(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');

atompos_next = 0;
% Reshape the atomic positions for Relax/MD
if S.RelaxFlag || S.MDFlag
	%atompos_next = transpose(S.lat_uvec) * S.Atoms';
	atompos_next = coordinateTransformation(S, S.Atoms, 'noncart2cart_dis');
	% update mapped Atomic position history, for keeping track of un-mapped
	% positions that are used in charge extrapolation
	S.Atoms_old = atompos_next;
	atompos_next = reshape(atompos_next',[],1);
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = initElectrondensity(S)
if (S.ForceCount == 1)
	% use sum of the atomic charge densities as rho_guess
	S.rho = S.rho_at;
else
	% perform charge extrapolation
	if ((S.RelaxFlag == 1 || S.MDFlag) && S.ForceCount >= 4)
		S.rho(:,1) = S.rho_at(:,1) + S.delta_rho_in_tp1;
		S.rho(S.rho(:,1) < 0,1) = S.xc_rhotol;
		% update spin up/down densities
		if S.nspin ~= 1
			rho_mag = S.rho(:,2) - S.rho(:,3);
			S.rho(:,2) = (S.rho(:,1) + rho_mag) * 0.5;
			S.rho(:,3) = (S.rho(:,1) - rho_mag) * 0.5;
		end
	end
	% need to scale the density
	S.rho = S.rho * abs(S.NegCharge/dot(S.W,S.rho(:,1)));
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = check_atomlocation(S)
% map atom positions back to the domain if atoms are outside the domain in
% periodic directions or throw an error in Dirichlet directions
if S.BCx == 1
	if (sum(S.Atoms(:,1) >= S.L1 | S.Atoms(:,1) < 0) > 0)
		error('Atom coordinates in the first lattice vector direction are out of cell');
	end
else
	S.Atoms(S.Atoms(:,1) >= S.L1,1) = S.Atoms(S.Atoms(:,1) >= S.L1,1) - S.L1;
	S.Atoms(S.Atoms(:,1) < 0,1) = S.Atoms(S.Atoms(:,1) < 0,1) + S.L1;
end

if S.BCy == 1
	if (sum(S.Atoms(:,2) >= S.L2 | S.Atoms(:,2) < 0) > 0)
		error('Atom coordinates in the second lattice vector direction are out of cell');
	end
else
	S.Atoms(S.Atoms(:,2) >= S.L2,2) = S.Atoms(S.Atoms(:,2) >= S.L2,2) - S.L2;
	S.Atoms(S.Atoms(:,2) < 0,2) = S.Atoms(S.Atoms(:,2) < 0,2) + S.L2;
end

if S.BCz == 1
	if (sum(S.Atoms(:,3) >= S.L3 | S.Atoms(:,3) < 0) > 0)
		error('Atom coordinates in the third lattice vector direction are out of cell');
	end
else
	S.Atoms(S.Atoms(:,3) >= S.L3,3) = S.Atoms(S.Atoms(:,3) >= S.L3,3) - S.L3;
	S.Atoms(S.Atoms(:,3) < 0,3) = S.Atoms(S.Atoms(:,3) < 0,3) + S.L3;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nndis = calculate_min_distance(S)
if S.n_atm == 1
    pairs = [1, 1];
else 
    pairs = nchoosek(1:S.n_atm,2);
end
p1 = S.Atoms(pairs(:,1),:);
p2 = S.Atoms(pairs(:,2),:);

dd = calculateDistance(p1(:,1),p1(:,2),p1(:,3),p2(:,1),p2(:,2),p2(:,3),S);
nndis = min(dd);
end







