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
mix_hist = S.MixingHistory; 
mix_p_pulay = S.PulayFrequency;
mix_param = S.MixingParameter; % beta
omega = S.MixingParameterSimple; % mixing parameter for simple mixing step

f_k = g_k - x_k;
if iter > 1
	f_km1 = S.mixing_hist_fkm1;
	x_km1 = S.mixing_hist_xkm1;
end

% store residual & iteration history
if iter > 1
	i_hist = mod(iter-2,mix_hist)+1;
	if (S.PulayRestartFreq > 0 && mod(iter, S.PulayRestartFreq) == 0)
		S.X = zeros(size(S.X)); S.F = zeros(size(S.F));
		S.X(:,1) = x_k - x_km1;
		S.F(:,1) = f_k - f_km1;
	else
		S.X(:,i_hist) = x_k - x_km1;
		S.F(:,i_hist) = f_k - f_km1;
	end
end

% apply Anderson extrapolation every p iters
if rem(iter,mix_p_pulay) == 0 && iter > 1
	% find weighted averages x_wavg, f_wavg
	[x_wavg, f_wavg] = andersonWtdAvg(x_k, f_k, S.X, S.F);
else 
	% simple mixing
	x_wavg = x_k; f_wavg = f_k;
end

% apply preconditioner to the weighted residual
if S.MixingPrecond > 0
	if S.nspin == 1
		% precondition total density
		if S.MixingPrecond == 1      % kerker precond
			S.Pf_guess = zeros(size(f_wavg)); % similar to using previous guess
			f_wavg = Kerker_Precond(S, f_wavg, S.precond_tol, 1000, S.Pf_guess);
			% S.Pf_guess = f_wavg; % update guess vector
		elseif S.MixingPrecond == 2  % resta precond
			S.Pf_guess = 1e-6*rand(length(f_wavg),length(S.precondcoeff_a));
			%S.precondcoeff_lambda_TF
			f_wavg = Resta_Precond(S, f_wavg, S.precond_tol, 1000, S.Pf_guess);
			%S.Pf_guess = repmat(f_wavg,1,length(S.precondcoeff_a)); % update guess vector
		elseif S.MixingPrecond == 3  % truncated kerker precond
			S.Pf_guess = 1e-6*rand(length(f_wavg),length(S.precondcoeff_a));
			%S.precondcoeff_lambda_TF
			f_wavg = TruncatedKerker_Precond(S, f_wavg, S.precond_tol, 1000, S.Pf_guess);
			%S.Pf_guess = repmat(f_wavg,1,length(S.precondcoeff_a)); % update guess vector
		end
	else
		f_Rho = f_wavg(1:S.N) + f_wavg(S.N+1:end);
		f_M = f_wavg(1:S.N) - f_wavg(S.N+1:end);
		% precondition total density
		if S.MixingPrecond == 1      % kerker precond
			Pf_Rho = Kerker_Precond(S, f_Rho, S.precond_tol, 1000, S.Pf_guess); 
			S.Pf_guess = Pf_Rho; % update guess vector
		elseif S.MixingPrecond == 2  % resta precond
			S.Pf_guess = 1e-6*rand(length(f_Rho),length(S.precondcoeff_a));
			Pf_Rho = Resta_Precond(S, f_Rho, S.precond_tol, 1000, S.Pf_guess); 
			%S.Pf_guess = repmat(Pf_Rho,1,length(S.precondcoeff_a)); % update guess vector
		elseif S.MixingPrecond == 3  % truncated kerker precond
			S.Pf_guess = 1e-6*rand(length(f_Rho),length(S.precondcoeff_a));
			Pf_Rho = TruncatedKerker_Precond(S, f_Rho, S.precond_tol, 1000, S.Pf_guess);
			%S.Pf_guess = repmat(Pf_Rho,1,length(S.precondcoeff_a)); % update guess vector
		end
		f_wavg = vertcat(Pf_Rho + f_M, Pf_Rho - f_M)/2;
	end
end    

% update new input iterate
if (rem(iter,mix_p_pulay) == 0)
	x_kp1 = x_wavg + mix_param * f_wavg;
else
	x_kp1 = x_wavg + omega * f_wavg;
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
function Pf = Kerker_Precond(S, f, tol, maxit, Pf_guess)
% KERKER_PRECOND applies kerker preconditioner.
%   KERKER_PRECOND effectively applies (Lap - k_TF^2)^-1 * Lap to f
%   by solving the linear system (Lap - k_TF^2) s = Lap*f. 

% lambda_sqr = lambda_TF * lambda_TF;
% 
% % these are fixed, so it should be possible to store them
% B = S.Lap_std - spdiags(lambda_sqr * ones(S.N,1),0,S.N,S.N);
% 
% Df = S.Lap_std * f;
% 
% Pf = aar(B,Df,Pf_guess,tol,maxit,[],[],[],[],[]);

Pf = RSfit_Precond(S,f,S.precondcoeff_a,S.precondcoeff_lambda_TF,...
				   S.precondcoeff_k,tol,maxit,Pf_guess);

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

