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
	[x_wavg, f_wavg] = andersonWtdAvg(x_k, f_k, S.X, S.F,S.nspden,S.MixingVariable);
else 
	% simple mixing
	x_wavg = x_k; f_wavg = f_k;
end

% calculate sum of all columns
f_wavg = reshape(f_wavg,[],S.nspden);
if S.spin_typ > 0
    sum_f = sum(f_wavg);
end

Pf = zeros(S.N,S.nspden);
% apply preconditioner to the residual of rho
if S.MixingPrecond > 0
    % precondition total density
    if S.MixingPrecond == 1      % kerker precond
        % S.Pf_guess = zeros(size(f_wavg)); % similar to using previous guess
        % precondition the residual of total density/potential
        k_TF = S.precond_kerker_kTF;
        idiemac = S.precond_kerker_thresh;
        Pf(:,1) = Kerker_Precond(S, f_wavg(:,1), amix, k_TF, idiemac, S.precond_tol, 1000, zeros(S.N,1));
    end
else
    Pf(:,1) = amix * f_wavg(:,1);     % mixing param is included in Pf
end

% apply preconditioner to the residual of magnetization
if S.MixingPrecondMag == 1
    k_TF_mag = S.precond_kerker_kTF_mag;
    idiemac_mag = S.precond_kerker_thresh_mag;
    for i = 2:S.nspden
        Pf(:,i) = Kerker_Precond(S, f_wavg(:,i), amix, k_TF_mag, idiemac_mag, S.precond_tol, 1000, zeros(S.N,1));
    end
else
    for i = 2:S.nspden
        Pf(:,i) = amix_mag * f_wavg(:,i);     % mixing param is included in Pf
    end
end

% this step makes scf converge slower when without spin
% don't know why
if S.spin_typ > 0
    % recover original sum
    sum_Pf = sum(Pf);
    shift = (sum_f - sum_Pf)/S.N;
    Pf = Pf + shift;
end

% get x_kp1
x_kp1 = x_wavg + reshape(Pf,[],1);

if S.MixingVariable == 0
	% due to inaccurate kerker solver, the density might have
	% slightly inaccuate integral, scale the density
    negrho_count = sum(x_kp1(1:S.N) < 0);
    if (negrho_count > 0)
        fprintf('\nDensity got negative\n\n');                
    end
    x_kp1(x_kp1(1:S.N) < 0) = 0;
    integral = S.W'*(x_kp1(1:S.N));
    x_kp1(1:S.N) = x_kp1(1:S.N) * (-S.NegCharge/integral);
end

% update the history vectors
S.mixing_hist_fkm1 = f_k;
S.mixing_hist_xkm1 = x_k;
S.mixing_hist_xk = x_kp1;
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

