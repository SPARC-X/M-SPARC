function [x_wavg, f_wavg] = andersonWtdAvg(x_k, f_k, X, F)
% @brief    ANDERSONWTDAVG finds the weighted averages
%           x_wavg := x_k - X*Gamma
%           f_wavg := f_k - F*Gamma
%           where Gamma = inv(F' * F) * F' * f_k, F' denotes the transpose 
%           of F.
%
% @param x_k   Current guess/input vector.
% @param f_k   Current preconditioned residual.
% @param X     Iterate histories, X = [x_{k-m+1}-x_{k-m}, ..., x_k-x_{k-1}].
% @param F     Residual histories, F = [f_{k-m+1}-f_{k-m}, ..., f_k-f_{k-1}].
%
% @authors  Qimen Xu <qimenxu@gatech.edu>
%           Abhiraj Sharma <asharma424@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @copyright (c) 2019 Material Physics & Mechanics Group, Georgia Tech

Gamma = pinv(F'*F)*(F'*f_k);  
x_wavg = x_k - X * Gamma;
f_wavg = f_k - F * Gamma;

end