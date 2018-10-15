%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program runs a sorted pairwise sorted QR deomposition
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q, R, p] = sorted_pairwiseQR(H, strt)

[nr, nc] = size(H);

R = zeros(nc);
Q = H;
p = 1 : 1 : nc;

sqnorm = zeros(1, nc/2);
for i = 1 : nc/2
    sqnorm(1, i) = Q(:, 2*i-1)' * Q(:, 2*i-1);
end

for i = 1 : nc/2
    [~, k] = min(sqnorm(1, i:nc/2));
    
    k = k + i - 1;

    R(:, [2*i-1, 2*k-1]) = R(:, [2*k-1, 2*i-1]);
    Q(:, [2*i-1, 2*k-1]) = Q(:, [2*k-1, 2*i-1]);
    p([2*i-1, 2*k-1]) = p([2*k-1, 2*i-1]);
    
    R(:, [2*i, 2*k]) = R(:, [2*k, 2*i]);
    Q(:, [2*i, 2*k]) = Q(:, [2*k, 2*i]);
    p([2*i, 2*k]) = p([2*k, 2*i]);

    sqnorm([i, k]) = sqnorm([k, i]);

    R(2*i-1, 2*i-1) = sqrt(sqnorm(i));
    Q(:, 2*i-1) = Q(:, 2*i-1) / R(2*i-1, 2*i-1);
    R(2*i, 2*i) = conj(R(2*i-1, 2*i-1));
   
    if strcmp(strt, 'Hbar')
        Q(1:(nr-nc)/2, 2*i) = -conj(Q((nr-nc)/2+1:(nr-nc), 2*i-1));
        Q((nr-nc)/2+1:(nr-nc), 2*i) = conj(Q(1:(nr-nc)/2, 2*i-1));
        for j = (nr - nc)/2 + 1 : nr/2
           Q(2*j-1, 2*i) = -conj(Q(2*j, 2*i-1));
           Q(2*j, 2*i) = conj(Q(2*j-1, 2*i-1));
        end
    elseif strcmp(strt, 'H')
       Q(1:nr/2, 2*i) = -conj(Q((nr/2+1):nr, 2*i-1));
       Q((nr/2+1):nr, 2*i) = conj(Q(1:nr/2, 2*i-1)); 
    end
       
    for n = i+1 : nc/2
        R(2*i-1,2*n-1) = Q(:, 2*i-1)' * Q(:, 2*n-1);
        R(2*i-1,2*n) = Q(:, 2*i-1)' * Q(:, 2*n);
        Q(:, 2*n-1) = Q(:, 2*n-1) - R(2*i-1,2*n-1)*Q(:,2*i-1) + conj(R(2*i-1,2*n)) * Q(:, 2*i);
        
        sqnorm(n) = sqnorm(n) - abs(R(2*i-1, 2*n-1))^2 - abs(R(2*i-1, 2*n))^2;

        R(2*i, 2*n-1) = -conj(R(2*i-1,2*n));
        R(2*i,2*n) = conj(R(2*i-1, 2*n-1));
       
        if strcmp(strt, 'Hbar')
            Q(1:(nr-nc)/2, 2*n) = -conj(Q((nr-nc)/2+1:(nr-nc), 2*n-1));
            Q((nr-nc)/2+1:(nr-nc), 2*n) = conj(Q(1:(nr-nc)/2, 2*n-1));
            for j = (nr - nc)/2 + 1 : nr/2
                Q(2*j-1, 2*n) = -conj(Q(2*j, 2*n-1));
                Q(2*j, 2*n) = conj(Q(2*j-1, 2*n-1));
            end 
        elseif strcmp(strt, 'H')
            Q(1:nr/2, 2*n) = -conj(Q((nr/2+1):nr, 2*n-1));
            Q((nr/2+1):nr, 2*n) = conj(Q(1:nr/2, 2*n-1));  
        end
    end
end

end