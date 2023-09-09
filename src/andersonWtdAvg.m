function [x_wavg, f_wavg] = andersonWtdAvg(x_k, f_k, X, F,nspden,opt)
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


FtF = compute_FtF(F,nspden,opt);
Ftf = compute_Ftf(F,f_k,nspden,opt);

Gamma = pinv(FtF)*Ftf;  
x_wavg = x_k - X * Gamma;
f_wavg = f_k - F * Gamma;

end

function FtF = compute_FtF(F,nspden,opt)
if nspden == 4
    m = size(F,2);
    FtF = zeros(m);
    for i = 1:m
        vec1 = F(:,i);
        for j = 1:m
            vec2 = F(:,j);
            FtF(i,j) = dotprod_nc(vec1,vec2,opt);
        end
    end
else
    FtF = F'*F;
end
end

function Ftf = compute_Ftf(F,f,nspden,opt)
if nspden == 4
    m = size(F,2);
    Ftf = zeros(m,1);
    for i = 1:m
        Ftf(i) = dotprod_nc(F(:,i),f,opt);
    end
else
    Ftf = F'*f;
end
end


function dotprod = dotprod_nc(vec1,vec2,opt)
if opt == 1 % potential
    N = length(vec1);
    p1 = (1:N/2);
    p2 = (1:N/2)+N/2;
    dotprod = sum(sum(vec1(p1).*vec2(p1)) + 2*sum(vec1(p2).*vec2(p2)));
else % density
    dotprod = 0.5*sum(vec1.*vec2);
end
end
