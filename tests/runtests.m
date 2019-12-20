% restoredefaultpath % warning: this will remove all the user-set path
addpath(genpath('../src'));

fnames = {'BC/H2O_sheet/H2O_Nx_28','BC/H2O_wire/H2O_wire','MeshConvergence/Si8-ONCV-0.7',...
		  'MeshConvergence/Si8-ONCV-0.6','MeshConvergence/Si8-ONCV-0.5','MeshConvergence/Si8-ONCV-0.4',...
		  'MeshConvergence/Si8-ONCV-0.3','MeshConvergence/Si8-ONCV-0.2', ...
		  'MD/NVE/Si8/Si8', 'MD/NVE/Si64/Si64' 'Nonorthogonal/Al/Al', ...
          'Nonorthogonal/Si/Si', 'Relax/Si/Si', 'Relax/SiH4/SiH4', 'Relax/vol_relax/Si/Si',...
          'Spin/Fe/Fe','Spin/O2/O2'};

for itest = 1:length(fnames)
    name = fnames{itest};
    clear S; clc; 
    % diary(strcat(name,'.log'));
    S = msparc(name);
    % save(strcat(name,'.mat'), 'S');
    % diary off;
end      


