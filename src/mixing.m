function [S,x_kp1] = mixing(S, g_k, x_k, iter)
% @brief    MIXING performs mixing of the previous SCF iterates to 
%           accelerate SCF convergence.
%
% @param g_k    Current output mixing iterate.
% @param x_in   Current input mixing iterate.
% @param iter   Iteration number.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech

[S,x_kp1] = Periodic_Pulay(S, g_k, x_k, iter);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S,x_kp1] = Periodic_Pulay(S, g_k, x_k, iter)
% periodic pulay mixing
% Ref: A Banerjee et al, Chemical Physics Letters 647 (2016) 31-35.
%      http://dx.doi.org/10.1016/j.cplett.2016.01.033
m = S.MixingHistory; 
p = S.PulayFrequency;
beta = S.MixingParameter; % beta
omega = S.MixingParameterSimple; % mixing parameter for simple mixing step
beta_mag = S.MixingParameterMag;
omega_mag = S.MixingParameterSimpleMag;

Pulay_mixing_flag = (rem(iter,p) == 0 && iter > 1);

if Pulay_mixing_flag   % paulay mixing
    amix = beta; amix_mag = beta_mag;
else                   % simple (linear) mixing
    amix = omega; amix_mag = omega_mag;
end

f_k = g_k - x_k;
if iter > 1
	f_km1 = S.mixing_hist_fkm1;
	x_km1 = S.mixing_hist_xkm1;
end

% store residual & iteration history
if iter > 1
	i_hist = mod(iter-2,m)+1;
	if (S.PulayRestartFlag ~= 0 && i_hist == 1)
		S.X = zeros(size(S.X)); S.F = zeros(size(S.F));
		S.X(:,1) = x_k - x_km1;
		S.F(:,1) = f_k - f_km1;
	else
		S.X(:,i_hist) = x_k - x_km1;
		S.F(:,i_hist) = f_k - f_km1;
	end
end

% apply Anderson extrapolation every p iters
if Pulay_mixing_flag
	% find weighted averages x_wavg, f_wavg
	[x_wavg, f_wavg] = andersonWtdAvg(x_k, f_k, S.X, S.F);
else 
	% simple mixing
	x_wavg = x_k; f_wavg = f_k;
end

sum_f_tot = 0.0;
sum_f_mag = 0.0;
if S.spin_typ ~= 0
    f_tot = f_wavg(1:S.N) + f_wavg(S.N+1:2*S.N);
    f_mag = f_wavg(1:S.N) - f_wavg(S.N+1:2*S.N);
    sum_f_tot = sum(f_tot);
    sum_f_mag = sum(f_mag);
else
    f_tot = f_wavg; % for spin-unpolarized calculations, f_tot is just f_wavg
    % f_mag is N/A for spin-unporlaized calculations
end

Pf = zeros(S.N,S.nspin);
% apply preconditioner to the weighted residual
if S.MixingPrecond > 0
    % precondition total density
    if S.MixingPrecond == 1      % kerker precond
        % S.Pf_guess = zeros(size(f_wavg)); % similar to using previous guess
        % precondition the residual of total density/potential
        k_TF = S.precond_kerker_kTF;
        idiemac = S.precond_kerker_thresh;
        Pf(:,1) = Kerker_Precond(S, f_tot, amix, k_TF, idiemac, S.precond_tol, 1000, zeros(S.N,1));
    end
else
    Pf(:,1) = amix * f_tot;     % mixing param is included in Pf
end

if S.spin_typ ~= 0
    if S.MixingPrecondMag ~= 0
        if S.MixingPrecondMag == 1
            k_TF_mag = S.precond_kerker_kTF_mag;
            idiemac_mag = S.precond_kerker_thresh_mag;
            Pf(:,2) = Kerker_Precond(S, f_mag, amix, k_TF_mag, idiemac_mag, S.precond_tol, 1000, zeros(S.N,1));
        end
    else
        Pf(:,2) = amix_mag * f_mag;     % mixing param is included in Pf
    end
end

if S.spin_typ ~= 0
    sum_Pf_tot = sum(Pf(:,1));
    shift_Pf_tot = (sum_f_tot - sum_Pf_tot)/S.N;
    Pf(:,1) = Pf(:,1) + shift_Pf_tot;
    
    sum_Pf_mag = sum(Pf(:,2));
    shift_Pf_mag = (sum_f_mag - sum_Pf_mag)/S.N;
    Pf(:,2) = Pf(:,2) + shift_Pf_mag;
end

% x_kp1 = x_wavg + mix_param * f_wavg;
if S.spin_typ == 0
    x_kp1 = x_wavg + Pf;
else
    x_kp1 = x_wavg + vertcat(Pf(:,1) + Pf(:,2), Pf(:,1) - Pf(:,2))/2;
end


if S.MixingVariable == 0
	% due to inaccurate kerker solver, the density might have
	% slightly inaccuate integral
	if S.MixingPrecond ~= 0
		% scale the density
		x_kp1 = x_kp1 * (-S.NegCharge/sum(S.W'*reshape(x_kp1,S.N,S.nspin),2));
	end
end

% update the history vectors
S.mixing_hist_fkm1 = f_k;
S.mixing_hist_xkm1 = x_k;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pf = Kerker_Precond(S, f, a, lambda_TF, idiemac, tol, maxit, Pf_guess)
% KERKER_PRECOND applies kerker preconditioner.
% For given function f, this function returns 
% Pf := a * (L - lambda_TF^2)^-1 * (L - idemac*lambda_TF^2)f, 
% where L is the discrete Laplacian operator, c is the 
% inverse of diemac (dielectric macroscopic constant).
% When c is 0, it's the original Kerker preconditioner.

% lambda_sqr = lambda_TF * lambda_TF;
% 
% % these are fixed, so it should be possible to store them
% B = S.Lap_std - spdiags(lambda_sqr * ones(S.N,1),0,S.N,S.N);
% 
% Df = S.Lap_std * f;
% 
% Pf = aar(B,Df,Pf_guess,tol,maxit,[],[],[],[],[]);

% Pf = RSfit_Precond(S,f,S.precondcoeff_a,S.precondcoeff_lambda_TF,...
% 				   S.precondcoeff_k,tol,maxit,Pf_guess);

Lf = S.Lap_std * f - (lambda_TF*lambda_TF*idiemac)*f;
B = S.Lap_std - spdiags(lambda_TF*lambda_TF * ones(S.N,1),0,S.N,S.N);
Pf = aar(B,Lf,Pf_guess,tol,maxit,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU);
if abs(lambda_TF) < 1E-14
    Pf = Pf - sum(Pf)/S.N;
end
Pf = a * Pf;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pf = Resta_Precond(S, f, tol, maxit, Pf_guess)
% RESTA_PRECOND applies real-space resta preconditioner.

Pf = RSfit_Precond(S,f,S.precondcoeff_a,S.precondcoeff_lambda_TF,...
				   S.precondcoeff_k,tol,maxit,Pf_guess);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pf = TruncatedKerker_Precond(S, f, tol, maxit, Pf_guess)
% TRUNCATEDKERKER_PRECOND applies real-space truncated-Kerker preconditioner.

Pf = RSfit_Precond(S,f,S.precondcoeff_a,S.precondcoeff_lambda_TF,...
				   S.precondcoeff_k,tol,maxit,Pf_guess);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pf = RSfit_Precond(S, f, a, lambda_sqr, k, tol, maxit, Pf_guess)
% RSFIT_PRECOND applies real-space preconditioner with any rational fit 
% coefficients.
%   RSFIT_PRECOND effectively applies sum_i (a_i*(Lap - k_TF2_i)^-1 * Lap + k*I) 
%   to f by solving the linear systems a_i*(Lap - k_TF2_i) s = Lap*f and 
%   summing the sols. To apply any preconditioner, simply perform a
%   rational curve fit to the preconditioner in fourier space and provide
%   the fit coeffs a(i), lambda_TF(i) and a constant k.

m = length(a); % number of terms

Df = S.Lap_std * f;

Pf = k * f;

isguess = size(Pf_guess,2) >= m; 
for i = 1:m
	B = S.Lap_std - spdiags(lambda_sqr(i) * ones(S.N,1),0,S.N,S.N);
	guess = []; if isguess, guess = Pf_guess(:,i); end
	Pf = Pf + a(i) * aar(B,Df,guess,tol,maxit,0.6,0.6,7,6,S.LapPreconL,S.LapPreconU);
	%[f1,~,~,~] = gmres(B,Df,50,tol,maxit,S.LapPreconL,S.LapPreconU,guess);
	%Pf = Pf + a(i) * f1;
end

Pf = real(Pf);
end

