function mass = atomdata_mass(element)
% @brief    Find atomic mass info. for given element type.
%
% @param element    The element name
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
% @ref     NIST (https:%physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=html)
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%

if (strcmp(element,'H'))
	mass = 1.007975; 
elseif (strcmp(element,'He'))
	mass = 4.002602; % avg min&max
elseif (strcmp(element,'Li'))
	mass = 6.9675; % avg min&max
elseif (strcmp(element,'Be'))
	mass = 9.0121831;
elseif (strcmp(element,'B'))
	mass = 10.8135; % avg min&max
elseif (strcmp(element,'C'))
	mass = 12.0106; % avg min&max
elseif (strcmp(element,'N'))
	mass = 14.006855; % avg min&max
elseif (strcmp(element,'O'))
	mass = 15.9994; % avg min&max
elseif (strcmp(element,'F'))
	mass = 18.998403163;
elseif (strcmp(element,'Ne'))
	mass = 20.1797;
elseif (strcmp(element,'Na'))
	mass = 22.98976928;
elseif (strcmp(element,'Mg'))
	mass = 24.3055; % avg min&max
elseif (strcmp(element,'Al'))
	mass = 26.9815385;
elseif (strcmp(element,'Si'))
	mass = 28.085; % avg min&max
elseif (strcmp(element,'P'))
	mass = 30.973761998;
elseif (strcmp(element,'S'))
	mass = 32.0675; % avg min&max
elseif (strcmp(element,'Cl'))
	mass = 35.4515; % avg min&max
elseif (strcmp(element,'Ar'))
	mass = 39.948;
elseif (strcmp(element,'K'))
	mass = 39.0983;
elseif (strcmp(element,'Ca'))
	mass = 40.078;
elseif (strcmp(element,'Sc'))
	mass = 44.955908;
elseif (strcmp(element,'Ti'))
	mass = 47.867;
elseif (strcmp(element,'V'))
	mass = 50.9415;
elseif (strcmp(element,'Cr'))
	mass = 51.9961;
elseif (strcmp(element,'Mn'))
	mass = 54.938044;
elseif (strcmp(element,'Fe'))
	mass = 55.845;
elseif (strcmp(element,'Co'))
	mass = 58.933194;
elseif (strcmp(element,'Ni'))
	mass = 58.6934;
elseif (strcmp(element,'Cu'))
	mass =  63.546; 
elseif (strcmp(element,'Zn'))
	mass =  65.38; 
elseif (strcmp(element,'Ga'))
	mass =  69.723; 
elseif (strcmp(element,'Ge'))
	mass = 72.630; 
elseif (strcmp(element,'As'))
	mass = 74.921595; 
elseif (strcmp(element,'Se'))
	mass = 78.971; 
elseif (strcmp(element,'Br'))
	mass = 79.904; % avg min&max
elseif (strcmp(element,'Kr'))
	mass = 83.798; 
elseif (strcmp(element,'Rb'))
	mass = 85.4678; 
elseif (strcmp(element,'Sr'))
	mass = 87.62; 
elseif (strcmp(element,'Y'))
	mass = 88.90584; 
elseif (strcmp(element,'Zr'))
	mass = 91.224; 
elseif (strcmp(element,'Nb'))
	mass = 92.90637; 
elseif (strcmp(element,'Mo'))
	mass = 95.95; 
elseif (strcmp(element,'Tc'))
	mass = 98; 
elseif (strcmp(element,'Ru'))
	mass = 101.07; 
elseif (strcmp(element,'Rh'))
	mass = 102.90550; 
elseif (strcmp(element,'Pd'))
	mass = 106.42; 
elseif (strcmp(element,'Ag'))
	mass = 107.8682; 
elseif (strcmp(element,'Cd'))
	mass = 112.414; 
elseif (strcmp(element,'In'))
	mass = 114.818; 
elseif (strcmp(element,'Sn'))
	mass = 118.710; 
elseif (strcmp(element,'Sb'))
	mass = 121.760; 
elseif (strcmp(element,'Te'))
	mass = 127.60; 
elseif (strcmp(element,'I'))
	mass = 126.90447; 
elseif (strcmp(element,'Xe'))
	mass = 131.293; 
elseif (strcmp(element,'Cs'))
	mass = 132.90545196; 
elseif (strcmp(element,'Ba'))
	mass = 137.327; 
elseif (strcmp(element,'La'))
	mass = 138.90547; 
elseif (strcmp(element,'Ce'))
	mass = 140.116; 
elseif (strcmp(element,'Pr'))
	mass = 140.90766; 
elseif (strcmp(element,'Nd'))
	mass = 144.242; 
elseif (strcmp(element,'Pm'))
	mass = 145; 
elseif (strcmp(element,'Sm'))
	mass = 150.36; 
elseif (strcmp(element,'Eu'))
	mass = 151.964; 
elseif (strcmp(element,'Gd'))
	mass = 157.25; 
elseif (strcmp(element,'Tb'))
	mass = 158.92535; 
elseif (strcmp(element,'Dy'))
	mass = 162.500; 
elseif (strcmp(element,'Ho'))
	mass = 164.93033; 
elseif (strcmp(element,'Er'))
	mass = 167.259; 
elseif (strcmp(element,'Tm'))
	mass = 168.93422; 
elseif (strcmp(element,'Yb'))
	mass = 173.054; 
elseif (strcmp(element,'Lu'))
	mass = 174.9668; 
elseif (strcmp(element,'Hf'))
	mass = 178.49; 
elseif (strcmp(element,'Ta'))
	mass = 180.94788; 
elseif (strcmp(element,'W'))
	mass = 183.84; 
elseif (strcmp(element,'Re'))
	mass = 186.207; 
elseif (strcmp(element,'Os'))
	mass = 190.23; 
elseif (strcmp(element,'Ir'))
	mass = 192.217; 
elseif (strcmp(element,'Pt'))
	mass = 195.084; 
elseif (strcmp(element,'Au'))
	mass = 196.966569; 
elseif (strcmp(element,'Hg'))
	mass = 200.592; 
elseif (strcmp(element,'Tl'))
	mass = 204.3835; % avg min&max
elseif (strcmp(element,'Pb'))
	mass = 207.2; 
elseif (strcmp(element,'Bi'))
	mass = 208.98040; 
elseif (strcmp(element,'Po'))
	mass = 209; 
elseif (strcmp(element,'At'))
	mass = 210; 
elseif (strcmp(element,'Rn'))
	mass = 222; 
elseif (strcmp(element,'Fr'))
	mass = 223; 
elseif (strcmp(element,'Ra'))
	mass = 226; 
elseif (strcmp(element,'Ac'))
	mass = 227; 
elseif (strcmp(element,'Th'))
	mass = 232.0377; 
elseif (strcmp(element,'Pa'))
	mass = 231.03588; 
elseif (strcmp(element,'U'))
	mass = 238.02891; 
elseif (strcmp(element,'Np'))
	mass = 237; 
elseif (strcmp(element,'Pu'))
	mass = 244; 
elseif (strcmp(element,'Am'))
	%mass = ; 
	fprintf('No atomic data for element %s, please provide in input file!\n', element);
elseif (strcmp(element,'Cm'))
	%mass = ; 
	fprintf('No atomic data for element %s, please provide in input file!\n', element);
elseif (strcmp(element,'Bk'))
	%mass = ; 
	fprintf('No atomic data for element %s, please provide in input file!\n', element);
elseif (strcmp(element,'Cf'))
	%mass = ; 
	fprintf('No atomic data for element %s, please provide in input file!\n', element);
elseif (strcmp(element,'Es'))
	%mass = ; 
	fprintf('No atomic data for element %s, please provide in input file!\n', element);
elseif (strcmp(element,'Fm'))
	%mass = ; 
	fprintf('No atomic data for element %s, please provide in input file!\n', element);
elseif (strcmp(element,'Md'))
	%mass = ; 
	fprintf('No atomic data for element %s, please provide in input file!\n', element);
elseif (strcmp(element,'No'))
	%mass = ; 
	fprintf('No atomic data for element %s, please provide in input file!\n', element);
elseif (strcmp(element,'Lw") == 0 || strcmp(element,"Lr'))
	%mass = ; 
	fprintf('No atomic data for element %s, please provide in input file!\n', element);
else
	fprintf('No atomic data for element %s, please provide in input file!\n', element);
end


end


