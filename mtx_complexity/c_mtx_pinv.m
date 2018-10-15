function c = c_mtx_pinv(M, N)

% when M > N, pinv(H) = inv(H' * H) * H'
% when N > M, pinv(H) = H' * inv(H * H')
if(M > N)
	c = c_mtx_mult(N, M, N) + c_mtx_inv(N)...
	 + c_mtx_mult(N, N, M);
else
	c = c_mtx_mult(M, N, M) + c_mtx_inv(M)...
	 + c_mtx_mult(N, M, M);
end