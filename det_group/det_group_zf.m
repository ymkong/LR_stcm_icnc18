%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong18apply] Applying Lattice Reduction Technique to Space-Time Coded Multiplexing Systems 
% This program implements group-zf detector
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mdl, m] = det_group_zf(act, mdl, m, lr_func, lr_name)

if(strcmp(act, 'updateH'))
    m.iterations = 0;
    m.basis_updates = 0;
    m.cmp_arithmetics = 0;
    m.arithmetics = 0;
    m.extra = 0;

elseif(strcmp(act, 'det'))
	Nr = mdl.Nr;
	Nt = mdl.Nt;
    H = mdl.chn_info_st.H;

    if(~(Nr == 4 && Nt == 4))
		disp('ERROR: group zf receiver undefined for current size.');
	else
		A = H(1:2, 1:2);
		B = H(1:2, 3:4); 
		C = H(3:4, 1:2);
		D = H(3:4, 3:4);

		W_zf = [inv(B), -inv(D); inv(A), -inv(C)];
		y_tilde = W_zf * m.y;

		s_hat = [];

		h_1 = inv(B) * A - inv(D) * C;
		Ht = [h_1(:, 1), h_1(:, 2); -conj(h_1(:, 2)), conj(h_1(:, 1))];
		y_1 = y_tilde(1:2, :);
		s_hat = [s_hat; Ht' * [y_1(:, 1); conj(y_1(:, 2))] / sum(diag(h_1' * h_1))];

		h_2 = inv(A) * B - inv(C) * D;
		Ht = [h_2(:, 1), h_2(:, 2); -conj(h_2(:, 2)), conj(h_2(:, 1))];
		y_2 = y_tilde(3:4, :);
		s_hat = [s_hat; Ht' * [y_2(:, 1); conj(y_2(:, 2))] / sum(diag(h_2' * h_2))];

		m.s_hat = s_hat;
        m.partb_cmp_arithmetics = c_mtx_inv(size(B, 1)) + c_mtx_inv(size(D, 1)) + c_mtx_inv(size(A, 1)) + c_mtx_inv(size(C, 1)) ...
        			+ c_mtx_mult(size(W_zf, 1), size(W_zf, 2), length(m.y)) ...
        			+ (c_mtx_mult(size(B, 1), size(A, 1), size(A, 2)) + c_mtx_mult(size(D, 1), size(C, 1), size(C, 2)) + c_mtx_add(size(B, 1), size(A, 2)) ...
        			+ c_mtx_mult(size(Ht, 2), size(Ht, 1), 1) + c_mtx_mult(size(h_1, 2), size(h_1, 1), size(h_1, 2) + size(h_1, 2) + size(Ht, 2))) * 2;
	end
end
end