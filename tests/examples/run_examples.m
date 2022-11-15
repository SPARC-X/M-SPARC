% This script launches four example tests
% Date: Nov 17,2019
% Author: Qimen Xu

% add src/ folder to path
% addpath(genpath('../../src'));
addpath(genpath(fullfile('..','src'))); % this works for all OSs

fnames = {'BaTiO3', 'H2O_sheet', 'H2O_wire', 'SiH4'};

t_s = tic;
for itest = 1:length(fnames)
    name = fnames{itest};
    clear S; clc; 
    diary(strcat(name,'.log'));
    %%%%%%%%%%%%%%%%%
    S = msparc(name);
    %%%%%%%%%%%%%%%%%
    % save(strcat(name,'.mat'), 'S');
    diary off;
end 
fprintf(' Total time for running all examples: %.3f s\n', toc(t_s)); 

