function c = c_mtx_qr(M, N)

% using the expression for SVD
% return number of complex operations
c = 4 * N^2 * M + 13 * M^3;