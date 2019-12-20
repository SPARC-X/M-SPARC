function S = chargeExtrapolation(S)
% @brief    Perform charge extrapolation to get better guess density for
%			MD/Relax.
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @ref: Dario Alfe, Computer Physics Communication 118 (1999) 31-33
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech

% update charge difference (rho_SC - rho_at) histories 
S.delta_rho_tm2 = S.delta_rho_tm1;
S.delta_rho_tm1 = S.delta_rho_t;

% find the difference between the current converged density and rho_at
% note that this is done before new rho_at is calculated
if (S.ForceCount >= 2)
	S.delta_rho_t = S.rho(:,1) - S.rho_at(:,1);
end

% find new extrapolation of charge difference for new rho_guess
% NOTE: this differs from delta_rho_t, delta_rho_tm1, delta_rho_tm2, since
% it is not actually a difference between self-consistent rho and rho_at.
% Instead, it's an extrapolation of the three
if (S.ForceCount >= 4)
	A = zeros(2,2); b = zeros(2,1);
	r_t_tm1   = S.atom_pos_t(:)   - S.atom_pos_tm1(:);
	r_tm1_tm2 = S.atom_pos_tm1(:) - S.atom_pos_tm2(:);
	r_tp1_t   = S.atom_pos_tp1(:) - S.atom_pos_t(:);
	A(1,1) = dot(r_t_tm1, r_t_tm1);
	A(1,2) = dot(r_t_tm1, r_tm1_tm2);
	A(2,2) = dot(r_tm1_tm2, r_tm1_tm2);
	A(2,1) = A(1,2);
	b(1)   = dot(r_tp1_t, r_t_tm1);
	b(2)   = dot(r_tp1_t, r_tm1_tm2);
	% find extrapolation coeff
	x = pinv(A) * b; % A might be singular
	alpha = x(1); beta = x(2);
	S.delta_rho_in_tp1 = (1+alpha) * S.delta_rho_t + ...
		(beta-alpha) * S.delta_rho_tm1 - beta * S.delta_rho_tm2;
	% S.delta_rho_in_tp1 = (1+alpha)*(S.dV_t/S.dV_tp1) * S.delta_rho_t + ...
	% 	(beta-alpha)*(S.dV_tm1/S.dV_tp1) * S.delta_rho_tm1 - beta*(S.dV_tm2/S.dV_tp1) * S.delta_rho_tm2;
end

% update atomic position histories 
S.atom_pos_tm2 = S.atom_pos_tm1;
S.atom_pos_tm1 = S.atom_pos_t;
S.atom_pos_t   = S.atom_pos_tp1;

% S.dV_tm2 = S.dV_tm1;
% S.dV_tm1 = S.dV_t;
% S.dV_t   = S.dV_tp1;

end

