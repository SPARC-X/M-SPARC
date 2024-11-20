function S = EvaluateEnergyAtom(S)
% @brief    Evaluates total energy of atom in the radial solver.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param S      Struct containing relevant fields
% @param iter   Iteration number
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================
Nd = S.Nd;
w = S.w;
int_scale = S.int_scale;

Eband = S.EigVal.*S.occ.matrix;
Eband = sum(Eband(:));
rho = S.rho(:,1);

Eelec_dc = -(S.r(2:Nd).^2).*(rho).*S.phi./int_scale(2:Nd);
Eelec_dc = 0.5*4*pi*w(2:Nd)*Eelec_dc;

if S.spinFlag == 0  % Spin Un-polarized
    if S.NLCC_flag
        Exc = (S.r(2:Nd).^2).*S.exc.*(rho+S.rho_Tilde)./int_scale(2:Nd);
    else
        Exc = (S.r(2:Nd).^2).*S.exc.*(rho)./int_scale(2:Nd);
    end
    
    Exc = 4*pi*w(2:Nd)*Exc;
    Exc_dc = (S.r(2:Nd).^2).*S.Vxc.*rho./int_scale(2:Nd);
    Exc_dc = 4*pi*w(2:Nd)*Exc_dc;
    
    % mGGA
    if S.ixc(3) == 1
        Eext_scan_dc =  (S.r(2:Nd).^2).*S.VxcScan3.*S.tau./int_scale(2:Nd);
        Eext_scan_dc = 4*pi*w(2:Nd)*Eext_scan_dc;
        Exc_dc = Exc_dc + Eext_scan_dc;
    end
else               % Spin Polarized
    if S.NLCC_flag
        Exc = (S.r(2:Nd).^2).*S.exc.*(rho+S.rho_Tilde)./int_scale(2:Nd);
    else
        Exc = (S.r(2:Nd).^2).*S.exc.*(rho)./int_scale(2:Nd);
    end
    
    Exc = 4*pi*w(2:Nd)*Exc;
    Exc_dc = sum(S.Vxc.*S.rho(:,2:3),2);
    Exc_dc = (S.r(2:Nd).^2).*Exc_dc./int_scale(2:Nd);
    Exc_dc = 4*pi*w(2:Nd)*Exc_dc;
    
    % mGGA
    if S.ixc(3) == 1
        Eext_scan_dc =  ([S.r(2:Nd) S.r(2:Nd)].^2).*S.VxcScan3.*S.tau(:, 2:3)./[int_scale(2:Nd) int_scale(2:Nd)];
        Eext_scan_dc = sum(4*pi*w(2:Nd)*Eext_scan_dc);
        Exc_dc = Exc_dc + Eext_scan_dc;
    end
end

if S.usefock < 2
    Etot = Eband + Exc - Exc_dc + Eelec_dc;
else
    Exc = Exc + S.Eex;
    Etot = Eband + Exc - Exc_dc + Eelec_dc - 2*S.Eex;
end

S.Eband = Eband;
S.Exc = Exc;
S.Exc_dc = Exc_dc;
S.Eelec_dc = Eelec_dc;
S.Etot = Etot;

fprintf("============================================================\n")
fprintf("\t\t   <strong>Energy (Hartree)</strong>\n")
fprintf("============================================================\n")
fprintf("Free energy\t\t\t: %0.10e\n",Etot)
fprintf("Band Structure energy\t\t: %0.10e\n",Eband)
fprintf("XC energy\t\t\t: %0.10e\n",Exc)
fprintf("XC Correction energy\t\t: % 0.10e\n",-Exc_dc)
fprintf("Hartree energy\t\t\t: %0.10e\n",Eelec_dc)

if S.usefock > 1
    fprintf("Exact exchange energy\t\t: %0.10e\n", S.Eex);
end

end

