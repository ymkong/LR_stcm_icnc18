%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong18apply] Applying Lattice Reduction Technique to Space-Time Coded Multiplexing Systems 
% This program implements group-mmse detector
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mdl, m] = det_group_mmse(act, mdl, m, lr_func, lr_name)

if(strcmp(act, 'updateH'))
    m.iterations = 0;
    m.basis_updates = 0;
    m.cmp_arithmetics = 0;
    m.arithmetics = 0;
    m.extra = 0;

elseif(strcmp(act, 'det'))
    Nr = mdl.Nr;
    Nt = mdl.Nt;
    Es = m.Es;
    sigma = m.sigma;
    H = mdl.chn_info_st.H; 
    
    partb_cmp_arithmetics = 0;

    s_hat = [];
    
    for k = 1 : Nt/2
        h_k = H(:, (2*k-1):2*k);
        W_k = h_k * h_k' * inv(H * H' + sigma^2/Es * eye(Nr));
        y_k = W_k * m.y;
        Ht = [h_k(:, 1), h_k(:, 2); -conj(h_k(:, 2)), conj(h_k(:, 1))];
        s_hat = [s_hat; Ht' * [y_k(:, 1); conj(y_k(:, 2))] / sum(diag(h_k' * h_k))];

        partb_cmp_arithmetics = partb_cmp_arithmetics + c_mtx_mult(size(h_k, 1), size(h_k, 2), size(h_k, 1)) + c_mtx_mult(size(h_k, 1), size(h_k, 1), size(H, 1)) + c_mtx_inv(size(H, 1))...
                 + c_mtx_mult(size(H, 1), size(H, 2), size(H, 1)) + size(H, 1)...
                 + c_mtx_mult(size(W_k, 1), size(W_k, 2), length(m.y))...
                 + c_mtx_mult(size(Ht, 2), size(Ht, 1), 1) + c_mtx_mult(size(h_k, 2), size(h_k, 1), size(h_k, 2) + size(h_k, 2) + size(Ht, 2));
    end

    m.s_hat = s_hat;
    m.partb_cmp_arithmetics = partb_cmp_arithmetics;
end