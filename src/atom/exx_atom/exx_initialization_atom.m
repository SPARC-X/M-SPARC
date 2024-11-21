function S = exx_initialization_atom(S)
% @brief    Initializes some parameters for exact exchange in atom.
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param S      Struct containing relevant fields
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================
if S.MAXIT_FOCK < 0
    S.MAXIT_FOCK = 20;
end
if S.MINIT_FOCK < 0
    S.MINIT_FOCK = 2;
end
if S.FOCK_TOL < 0
    S.FOCK_TOL = 0.2* S.SCF_tol;
end
if S.SCF_tol_init < 0
    S.SCF_tol_init = max(10*S.FOCK_TOL,1e-3);
end

if S.xc == 41
    S.hyb_mixing = 0.25;
end
end

