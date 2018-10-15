%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program runs a pairwise QR deomposition without sorting
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q, R] = pairwiseQR(H, strt)

[nr, nc] = size(H);

R = zeros(nc);
Q = H;

for i = 1 : nc/2
   R(2*i-1, 2*i-1) = norm(Q(:, 2*i-1));
   Q(:, 2*i-1) = Q(:, 2*i-1) / R(2*i-1, 2*i-1);
   R(2*i, 2*i) = R(2*i-1, 2*i-1);
   
   if strcmp(strt, 'C')
       for j = 1 : nr/2
          Q(2*j-1, 2*i) = -conj(Q(2*j, 2*i-1));
          Q(2*j, 2*i) = conj(Q(2*j-1, 2*i-1));
       end
   elseif strcmp(strt, 'Hbar')
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
       R(2*i, 2*n-1) = -conj(R(2*i-1,2*n));
       R(2*i,2*n) = conj(R(2*i-1, 2*n-1));
       
       if strcmp(strt, 'C')
           for j = 1 : nr/2
              Q(2*j-1, 2*n) = -conj(Q(2*j, 2*n-1));
              Q(2*j, 2*n) = conj(Q(2*j-1, 2*n-1));
           end
       elseif strcmp(strt, 'Hbar')
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