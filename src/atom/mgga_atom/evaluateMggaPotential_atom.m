function [VmGGA] = evaluateMggaPotential_atom(S, l, x, spin)
% @brief    Calculates the meta GGA potential for each "l" channel and spin.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param S      Struct containing relevant fields
% @param iter   Iteration number
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================
% spin =  0.5 for up
% spin = -0.5 for down

if S.spinFlag == 0
    V3 = S.VxcScan3;
else
    if spin == 0.5
        V3 = S.VxcScan3(:,1);
    elseif spin == -0.5
        V3 = S.VxcScan3(:,2);
    end
end

Nd = S.Nd;
r = S.r(2:Nd);
D = S.Gradient.matrix;
D = D(2:Nd,2:Nd);

%% Eigs way

if l~=0
    term1 = (l*(l+1)*(V3./(r.^2))).*x;
else
    term1 = zeros(size(x));
end

term2 = - (1./(r)).*(D*(r.*V3.*D*x) - D*(V3.*x));
VmGGA = 0.5*(term1 + term2);

end

