%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program runs a sorted QR deomposition
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q, R, p] = sorted_QR(H)

[Nr, Nt] = size(H);
R = zeros(Nt, Nt);
Q = H;
p = 1 : Nt;

for i = 1 : Nt
    [val, ind] = min(diag(Q(:, i : Nt)' * Q(:, i : Nt)));
    ind = ind + i - 1;
    Q(:, [i, ind]) = Q(:, [ind, i]);
    R(:, [i, ind]) = R(:, [ind, i]);
    p(:, [i, ind]) = p(:, [ind, i]);
    
    R(i, i) = sqrt(Q(:, i)' * Q(:, i));
    Q(:, i) = Q(:, i) ./ R(i, i);
    
    for j = (i + 1) : Nt
        R(i, j) = Q(:, i)' * Q(:, j);
        Q(:, j) = Q(:, j) - R(i, j) * Q(:, i);
    end
    
end

end