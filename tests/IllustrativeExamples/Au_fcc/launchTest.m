% add source/ folder to path
addpath(genpath('../../src'));

% lauch tests

%lat = 7.50 : 0.05 : 7.80;
lat = 7.00 : 0.05 : 8.00;
%lat = [7.00 : 0.05 : 7.45, 7.85 : 0.05: 8.00];
fnames = cell(1,length(lat)); % empty cell
for i = 1:length(lat)
	fnames{i} = num2str(sprintf('lat%.2f/Au_fcc',lat(i)));
end

t_s = tic;
for itest = 1:length(fnames)
	name = fnames{itest};
	clear S; clc; 
	diary(strcat(name,'.log'));
	%%%%%%%%%%%%%%%%%
	S = msparc(name,1); % 1 is to turn on parfor over k-points
	%%%%%%%%%%%%%%%%%
	% save data
	occ = S.occ;
	EigVal = S.EigVal;
	rho = S.rho;
	Efermi = S.lambda_f;
	save('Au_fcc.mat','occ','EigVal','rho','Efermi');
	diary off;
end 
fprintf(' Total time for running all examples: %.3f s\n', toc(t_s)); 


