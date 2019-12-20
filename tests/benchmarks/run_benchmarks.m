% restoredefaultpath % warning: this will remove all the user-set path
addpath(genpath('../../src'));

fnames = {'Al_Vacancy/Al', 'Au_13/Au_13', 'Au_FCC/Au', 'BaTiO3/BaTiO3', ...
          'C_511/C', 'CO/CO', 'Cu_FCC/Cu_fcc', 'Fe_411/Fe', 'H2O_sheet/H2O_sheet',...
          'H2O_wire/H2O_wire','MoS2/MoS2', 'Na_BCC/Na_bcc', 'O2/O2', 'SiH4/SiH4'};

for itest = 1:length(fnames)
    name = fnames{itest};
    clear S; clc; diary(strcat(name,'.log'));
    S = msparc(name);
    % save(strcat(name,'.mat'), 'S');
    diary off;
end 

