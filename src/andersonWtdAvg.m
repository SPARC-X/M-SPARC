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
            FtF(i,j) = dotprodm_v(vec1,vec2,nspden,opt);
        end
    end
else
    FtF = F'*F;
end
end

function Ftf = compute_Ftf(F,f,nspden,opt)
m = size(F,2);
Ftf = zeros(m,1);
for i = 1:m
    Ftf(i) = dotprodm_v(F(:,i),f,nspden,opt);
end
end


function dotprod = dotprodm_v(vec1,vec2,nspden,opt)
vec1 = reshape(vec1,[],nspden);
vec2 = reshape(vec2,[],nspden);
if nspden == 4
    if opt == 1 % potential
        dotprod = sum(sum(vec1(:,1:2).*vec2(:,1:2)) + 2*sum(vec1(:,3:4).*vec2(:,3:4)));
    else % density
        dotprod = 0.5*sum(vec1(:).*vec2(:));
    end
else
    dotprod = sum(vec1(:).*vec2(:));
end
end
