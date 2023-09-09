function S = electronDensity(S)
% @brief    Calculate new density based on the new states.
%
% @param S  A struct that contains the relevant fields.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Xin Jing <xjing30@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech
%

rho_d = zeros(S.N,S.nspinor);
for spinor = 1:S.nspinor
    ndrange = (1+(spinor-1)*S.N:spinor*S.N);
    nsshift = (spinor-1)*S.tnkpt*(S.spin_typ == 1);
    for kpt =1:S.tnkpt
         rho_d(:,spinor) = rho_d(:,spinor) + S.occfac*S.wkpt(kpt)* sum( S.psi(ndrange,:,kpt).*conj(S.psi(ndrange,:,kpt)).*S.occ(:,kpt+nsshift)',2);
    end
end
rho_d = real(rho_d);
rhotot = sum(rho_d,2);

if S.spin_typ == 0
    S.rho = rhotot;
elseif S.spin_typ == 1
    S.rho = [rhotot rho_d];
    S.mag = rho_d(:,1) - rho_d(:,2);
elseif S.spin_typ == 2
    spinor1 = (1:S.N);
    spinor2 = (S.N+1:2*S.N);
    rho_od = zeros(S.N,1);
    for kpt =1:S.tnkpt
         rho_od = rho_od + S.occfac*S.wkpt(kpt)* sum( S.psi(spinor1,:,kpt).*conj(S.psi(spinor2,:,kpt)).*S.occ(:,kpt)',2);
    end
    mx = 2*real(rho_od);
    my = -2*imag(rho_od);
    mz = rho_d(:,1) - rho_d(:,2);
    magnorm = sqrt(mx.^2 + my.^2 + mz.^2);
    S.rho = [rhotot 0.5*(rhotot+magnorm)  0.5*(rhotot-magnorm)];
    S.mag = [magnorm mx my mz];
end

end


