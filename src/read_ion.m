function S = read_ion(S, filename)
% @brief    READ_ION reads the filename.ion file.
%
% @param fname  The ion filename, without suffix.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%

fprintf(' Reading .ion file ...\n');

% open .ion file
fid1=fopen(strcat(filename,'.ion'),'r');
if (fid1 == -1)
    error('\n Cannot open file "%s.ion"\n',filename);
end

% first identify total number of atom types
typcnt = 0;
while (~feof(fid1))
    %fscanf(fid1,"%s",str);
    C_inpt = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
    str = char(C_inpt{:});
    if (strcmp(str, 'ATOM_TYPE:'))
        typcnt = typcnt + 1;
    end
    % skip current line
    textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);
end

fprintf(' Number of atom types : %d\n',typcnt);

if (typcnt < 1)
    error('\n Please provide at least one type of atoms!\n');
end

S.n_typ = typcnt;

%--------------------------------------------------------------------------

% reset file pointer to the start of the file
fseek(fid1, 0, 'bof'); % equivalent ot frewind(fid1)

% find total number of atoms
n_atom = 0;    % totoal num of atoms
typcnt = 0;    % atom type count
atmcnt_cum = zeros(S.n_typ+1,1);

S.Atm = repmat(struct([]), S.n_typ, 1);
while (~feof(fid1))
    C_inpt = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
    str = char(C_inpt{:}); % for R2016b and later, can use string()
    if (strcmp(str, 'ATOM_TYPE:'))
        typcnt = typcnt + 1;
        C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
        S.Atm(typcnt).typ = char(C_param{:});
    elseif (strcmp(str, 'N_TYPE_ATOM:'))
        C_inpt = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
        S.Atm(typcnt).n_atm_typ = C_inpt{1};
        S.Atm(typcnt).coords = zeros(S.Atm(typcnt).n_atm_typ, 3);
        S.Atm(typcnt).psdfname = 'undefined';
        S.Atm(typcnt).Mass   = 1e9; % set the default later
        S.Atm(typcnt).lloc   = 4; % default is 4 for oncv
        S.Atm(typcnt).psptyp = 1; % default is psp8 format
        S.Atm(typcnt).mag    = zeros(S.Atm(typcnt).n_atm_typ,3);
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);
        n_atom = n_atom + S.Atm(typcnt).n_atm_typ ;
        atmcnt_cum(typcnt+1) = n_atom;
    else
        % skip current line
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);
    end
end

fprintf(' Total number of atoms: %d\n',n_atom);

if (n_atom < 1)
    error('\n Please provide at least one atom!\n');
end

S.n_atm = n_atom;

S.Atoms = zeros(n_atom,3);

%--------------------------------------------------------------------------

% set default atom relax constraints to be all on (move in all DOFs)
S.mvAtmConstraint = ones(n_atom,3);

% set default atomic masses based on atom types (for MD)
if (1)  % atomic mass is only needed for MD
    for ityp = 1:S.n_typ
        % first identify element type
        elemType = S.Atm(ityp).typ;
        % remove the number appending the element type name if any
        elemType = elemType(isstrprop(elemType, 'alpha')); % e.g., Si2 -> Si
        % S.Atm(ityp).element = elemType;
        % find default atomic mass
        S.Atm(ityp).Mass = atomdata_mass(elemType);
        fprintf(' Default atomic mass for %s is %f\n',elemType,S.Atm(ityp).Mass);
    end
end

%--------------------------------------------------------------------------

% reset temp var
typcnt = 0;
atmcnt_coord = 0;

% allocate the size of the Isfrac vector which stores the coordinate type of each atom type
S.IsFrac = zeros(S.n_typ,1);
S.IsSpin = zeros(S.n_typ,1);

% reset file pointer to the start of the file
fseek(fid1, 0, 'bof'); % equivalent ot frewind(fid1)

while (~feof(fid1))
    C_inpt = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
    str = char(C_inpt{:}); % for R2016b and later, can use string()

    % enable commenting with '#'
    if (isempty(str) || str(1) == '#')
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
        continue;
    end

    if (strcmp(str, 'ATOM_TYPE:'))
        typcnt = typcnt + 1;
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
    elseif (strcmp(str, 'N_TYPE_ATOM:'))
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
    elseif (strcmp(str, 'COORD:'))
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
        %fprintf(' typcnt = %d, \n',typcnt);
        typ_atm_count = 0; % atom count for this type
        S.Atm(typcnt).coords = zeros(S.Atm(typcnt).n_atm_typ,3);
        for i = 1:S.Atm(typcnt).n_atm_typ
            typ_atm_count = typ_atm_count + 1;
            atmcnt_coord = atmcnt_coord + 1;
            C_param = textscan(fid1,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
            S.Atm(typcnt).coords(typ_atm_count,1) = C_param{1};
            S.Atm(typcnt).coords(typ_atm_count,2) = C_param{2};
            S.Atm(typcnt).coords(typ_atm_count,3) = C_param{3};
            S.Atoms(atmcnt_coord,:) = S.Atm(typcnt).coords(typ_atm_count,:);
            textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
        end
    elseif (strcmp(str, 'COORD_FRAC:'))
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
        %fprintf(' typcnt = %d, \n',typcnt);
        typ_atm_count = 0; % atom count for this type
        for i = 1:S.Atm(typcnt).n_atm_typ
            typ_atm_count = typ_atm_count + 1;
            atmcnt_coord = atmcnt_coord + 1;
            C_param = textscan(fid1,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
            S.Atm(typcnt).coords(typ_atm_count,1) = C_param{1} * S.L1;
            S.Atm(typcnt).coords(typ_atm_count,2) = C_param{2} * S.L2;
            S.Atm(typcnt).coords(typ_atm_count,3) = C_param{3} * S.L3;
            S.Atoms(atmcnt_coord,:) = S.Atm(typcnt).coords(typ_atm_count,:);
            textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
        end
        S.IsFrac(typcnt) = 1;
    elseif (strcmp(str, 'SPIN:'))
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
        %fprintf(' typcnt = %d, \n',typcnt);
        typ_atm_count = 0; % atom count for this type
        for i = 1:S.Atm(typcnt).n_atm_typ
            typ_atm_count = typ_atm_count + 1;
            C_param = textscan(fid1,'%f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
            num_non_nan = sum(cellfun(@(x) sum(~isnan(x)), C_param));
            if num_non_nan == 1
                S.Atm(typcnt).mag(typ_atm_count,1) = 0;
                S.Atm(typcnt).mag(typ_atm_count,2) = 0;
                S.Atm(typcnt).mag(typ_atm_count,3) = C_param{1};
            elseif num_non_nan == 3
                S.Atm(typcnt).mag(typ_atm_count,1) = C_param{1};
                S.Atm(typcnt).mag(typ_atm_count,2) = C_param{2};
                S.Atm(typcnt).mag(typ_atm_count,3) = C_param{3};
            else
                error("Please specify either spin in Z direction or spin in x, y, z directions.");
            end
            textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
        end
        S.IsSpin(typcnt) = 1;
    elseif (strcmp(str, 'RELAX:'))
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
        atmcnt_relax = atmcnt_cum(typcnt);
        %typ_atm_count = 0; % atom count for this type
        for i = 1:S.Atm(typcnt).n_atm_typ
            % typ_atm_count = typ_atm_count + 1;
            atmcnt_relax  = atmcnt_relax + 1;
            C_param = textscan(fid1,'%d %d %d',1,'delimiter',' ','MultipleDelimsAsOne',1);
            S.mvAtmConstraint(atmcnt_relax,1) = C_param{1,1};
            S.mvAtmConstraint(atmcnt_relax,2) = C_param{1,2};
            S.mvAtmConstraint(atmcnt_relax,3) = C_param{1,3};
            textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
        end
        fprintf('WARNING: All atoms will be forced to relax\n');
    elseif (strcmp(str, 'PSEUDO_POT:'))
        C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
        S.Atm(typcnt).psdfname = char(C_param{:});
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
        S.is_default_psd = 0; % switch off default psedopots
        fprintf(' pseudo_dir # %d = %s\n',typcnt,S.Atm(typcnt).psdfname);
    elseif (strcmp(str, 'ATOMIC_MASS:'))
        C_param = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
        S.Atm(typcnt).Mass = C_param{1};
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
    elseif ( isnumeric(str(1)) )
        error(['\nPlease specify the identifier before numbers!\n', ...
            'Reminder: check if the number of atoms specified is inconsistent\n' ...
            '          with the number of coordinates provided\n']);
    elseif (strcmp(str, 'HUBBARD:')) % Hubbard
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
    elseif (strcmp(str, 'U_ATOM_TYPE:')) % Hubbard
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
    elseif (strcmp(str, 'CORE_FLAG:')) % Hubbard
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
    elseif (strcmp(str, 'U_VAL:')) % Hubbard
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
    else
        fprintf('\nCannot recognize input flag: "%s"',str);
    end
end

if (atmcnt_coord ~= n_atom)
    error(['The number of coordinates provided is inconsistent '...
        'with the given number of atoms!']);
end

fprintf('\n COORD:\n');
disp(S.Atoms);

%--------------------------------------------------------------------------

%fprintf(' RELAX:\n');
%disp(S.mvAtmConstraint);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Hubbard Part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Sayan Bhowmik (sbhowmik9@gatech.edu)

S.hubbard_flag = 0;
% reset file pointer to the start of the file
fseek(fid1, 0, 'bof'); % equivalent ot frewind(fid1)

U_typ_flags = zeros(S.n_typ, 1);
typcnt_U = 0;
n_atom_U = 0;

while (~feof(fid1))
    C_inpt = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
    str = char(C_inpt{:}); % for R2016b and later, can use string()

    % enable commenting with '#'
    if (isempty(str) || str(1) == '#')
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
        continue;
    end

    if (strcmp(str, 'HUBBARD:'))
        disp('Found Hubbard Block');
        S.hubbard_flag = 1;
        textscan(fid1,'%s',1,'delimiter','\n','MultipleDelimsAsOne',0);  % skip current line
    elseif (strcmp(str, 'U_ATOM_TYPE:'))
        C_param = textscan(fid1,'%s',1,'delimiter',' ','MultipleDelimsAsOne',1);
        typ = char(C_param{:});
        id_OG = -1;

        % Extract info from original struct
        for J_a = 1 : S.n_typ
            if (strcmp(typ,S.Atm(J_a).typ))
                id_OG = J_a;
            end
        end

        if id_OG == -1
            error('This atom type does not exist in the original types, cannot apply U.');
        else
            typcnt_U = typcnt_U + 1;
        end

        U_typ_flags(id_OG) = 1;
        S.AtmU(typcnt_U).n_atm_typ = S.Atm(id_OG).n_atm_typ;
        n_atom_U = n_atom_U + S.AtmU(typcnt_U).n_atm_typ;
        S.AtmU(typcnt_U).coords = S.Atm(id_OG).coords;
        S.IsFracU(typcnt_U) = S.IsFrac(id_OG);
        S.AtmU(typcnt_U).hub_core = 0;
    elseif (strcmp(str, 'CORE_FLAG:'))
        C_inpt = textscan(fid1,'%f',1,'delimiter',' ','MultipleDelimsAsOne',1);
        if C_inpt{1} == 1
            S.AtmU(typcnt_U).hub_core = 1;
        end
    elseif (strcmp(str, 'U_VAL:'))
        C_inpt = textscan(fid1,'%f %f %f %f',1,'delimiter',' ','MultipleDelimsAsOne',1);
        S.AtmU(typcnt_U).Uval(1) = C_inpt{1};
        S.AtmU(typcnt_U).Uval(2) = C_inpt{2};
        S.AtmU(typcnt_U).Uval(3) = C_inpt{3};
        S.AtmU(typcnt_U).Uval(4) = C_inpt{4};

        for u_i = 1:4
            % Check if the current entry is NaN
            if isnan(S.AtmU(typcnt_U).Uval(u_i))
                % If it is NaN, throw error
                error("Please input U values in the correct format - 4 values corresponding to s p d f orbitals.")
            end
        end
    end
end

S.n_typ_U = typcnt_U;
if S.n_typ_U == 0 && S.hubbard_flag == 1
    error("Please specify atoms where you want to apply hubbard corrections.")
elseif S.n_typ_U > 0 
    S.hubbard_flag = 1;
    disp('Found Hubbard Block');
end

S.U_typ_flags = U_typ_flags;

if S.hubbard_flag == 1
    fprintf(' Total number of atoms having U corrections: %d\n',n_atom_U);

    S.n_atm_U = n_atom_U;

    S.AtomsU = zeros(n_atom_U,3);

    atmcnt_coord_U = 0;

    for U_typ = 1 : S.n_typ_U
        for J_a = 1 : S.AtmU(U_typ).n_atm_typ
            atmcnt_coord_U = atmcnt_coord_U + 1;
            S.AtomsU(atmcnt_coord_U,:) = S.AtmU(U_typ).coords(J_a,:);
        end
    end

    if (atmcnt_coord_U ~= n_atom_U)
    	error(['The number of coordinates provided is inconsistent '...
            'with the given number of atoms needing U corrections!']);
    end

    fprintf('\n U-corrections needed for following atom coordinates:\n COORD:\n');
    disp(S.AtomsU);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Close file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(fid1);

% find total mass
S.TotalMass = 0;
for ityp = 1:S.n_typ
    S.TotalMass = S.TotalMass + S.Atm(ityp).Mass * S.Atm(ityp).n_atm_typ;
end
end