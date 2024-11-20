function S = mixingHubbard(S, iter)
% @brief    Performs mixing of the occupation matrices for each atom.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param S      Struct containing relevant fields
% @param iter   Iteration number
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================

% loop over all atoms having U correction
for J = 1 : S.n_atm_U
    angnum = S.AtomU(J).angnum;
    rho_mn_out = zeros(angnum*angnum*S.nspinor,1);
    for spinor = 1 : S.nspinor
        shift = (spinor - 1)*angnum*angnum;
        copy_rho_mn = S.AtomU(J).rho_mn(:,:,spinor);
        rho_mn_out(1 + shift : angnum*angnum + shift) = copy_rho_mn(:);
    end
    
    % mix occupation matrix
    [S, rho_mn_out] = rho_mn_mix(S, rho_mn_out, S.AtomU(J).mixing_hist_xk, iter, J);

    counter = 1;
    for spinor = 1 : S.nspinor
        for m2 = 1 : angnum
            for m1 = 1 : angnum
                S.AtomU(J).rho_mn(m1,m2,spinor) = rho_mn_out(counter);
                counter = counter + 1;
            end
        end
        S.AtomU(J).rho_mn(:,:,spinor) = 0.5* (S.AtomU(J).rho_mn(:,:,spinor) ...
            + S.AtomU(J).rho_mn(:,:,spinor)');
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S, x_kp1] = rho_mn_mix(S, g_k, x_k, iter, atm_id)
% periodic pulay mixing of the occupation matrices of Atom atm_id using the
% same Gamma from density mixing performed
m = S.MixingHistory; 
p = S.PulayFrequency;
beta = S.MixingParameter; % beta
omega = S.MixingParameterSimple; % mixing parameter for simple mixing step

Pulay_mixing_flag = (rem(iter,p) == 0 && iter > 1);

if Pulay_mixing_flag   % paulay mixing
    amix = beta; 
else                   % simple (linear) mixing
    amix = omega;
end

f_k = g_k - x_k;
if iter > 1
	f_km1 = S.AtomU(atm_id).mixing_hist_fkm1;
	x_km1 = S.AtomU(atm_id).mixing_hist_xkm1;
end

% store residual & iteration history
if iter > 1
	i_hist = mod(iter-2,m)+1;
	if (S.PulayRestartFlag ~= 0 && i_hist == 1)
		S.AtomU(atm_id).X = zeros(size(S.AtomU(atm_id).X)); 
        S.AtomU(atm_id).F = zeros(size(S.AtomU(atm_id).F));
		S.AtomU(atm_id).X(:,1) = x_k - x_km1;
		S.AtomU(atm_id).F(:,1) = f_k - f_km1;
	else
		S.AtomU(atm_id).X(:,i_hist) = x_k - x_km1;
		S.AtomU(atm_id).F(:,i_hist) = f_k - f_km1;
	end
end

% apply Anderson extrapolation every p iters using the Gamma from density
% mixing
if Pulay_mixing_flag
    f_wavg = f_k - S.AtomU(atm_id).F * S.mix_Gamma;
    x_wavg = x_k - S.AtomU(atm_id).X * S.mix_Gamma;
else
    f_wavg = f_k; x_wavg = x_k;
end

Pf = amix * f_wavg;

% Get x_kp1
x_kp1 = x_wavg + Pf;

% update the history vectors
S.AtomU(atm_id).mixing_hist_fkm1 = f_k;
S.AtomU(atm_id).mixing_hist_xkm1 = x_k;
S.AtomU(atm_id).mixing_hist_xk = x_kp1;
end