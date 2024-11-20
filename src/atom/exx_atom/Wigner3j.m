function W3j = Wigner3j(j1,j2,j)
% @ brief     This function returns values for the special
%             case of Wigner 3j symbol [j1 j2 j;0 0 0]
%
% @authors  Sayan Bhowmik <sbhowmik9@gatech.edu>
%           Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
%
% @param S      Struct containing relevant fields
%
% @copyright (c) 2024 Material Physics & Mechanics Group, Georgia Tech
%===============================================================================
J = j1 + j2 + j;

if mod(J,2) == 1
    W3j = 0;
else
    g = J/2;
    W3j = (-1)^g;
    W3j = W3j*sqrt(factorial(2*g-2*j1)*factorial(2*g-2*j2)*factorial(2*g-2*j)/factorial(2*g+1));
    W3j = W3j*factorial(g)/(factorial(g-j1)*factorial(g-j2)*factorial(g-j));
end
end

