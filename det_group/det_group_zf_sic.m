%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong18apply] Applying Lattice Reduction Technique to Space-Time Coded Multiplexing Systems 
% This program implements group-zf-sic detector
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mdl, m] = det_group_zf_sic(act, mdl, m, lr_func, lr_name)

if(strcmp(act, 'updateH'))
    m.iterations = 0;
    m.basis_updates = 0;
    m.cmp_arithmetics = 0;
    m.arithmetics = 0;
    m.extra = 0;

elseif(strcmp(act, 'det'))
    Nt = mdl.Nt;
    H = mdl.chn_info_st.H; 
    
    partb_cmp_arithmetics = 0;

    [Q, R] = qr(H, 0);
    y_tilde = Q' * m.y;

    partb_cmp_arithmetics = partb_cmp_arithmetics + c_mtx_qr(size(H, 1), size(H, 2)) + c_mtx_mult(size(Q, 2), size(Q, 1), length(m.y));

    s_hat = zeros(Nt, 1);
    for k = Nt / 2 : -1 : 1
        Rt = R((2*k-1):2*k, (2*k-1):2*k);
        yt = y_tilde((2*k-1):2*k, :);

        Ht = [Rt(:, 1), Rt(:, 2); -conj(Rt(:, 2)), conj(Rt(:, 1))];
        s_hat((2*k-1):2*k) = Ht' * [yt(:, 1); conj(yt(:, 2))] / sum(diag(Rt' * Rt));

        partb_cmp_arithmetics = partb_cmp_arithmetics + c_mtx_mult(size(Ht, 2), size(Ht, 1), 1) + c_mtx_mult(size(Rt, 2), size(Rt, 1), size(Rt, 2) + size(Rt, 2) + size(Ht, 2));

        y_tilde = y_tilde - R(:, (2*k-1):2*k) * [s_hat((2*k-1):2*k), [conj(s_hat(2*k)); -conj(s_hat(2*k-1))]];

        partb_cmp_arithmetics = partb_cmp_arithmetics + c_mtx_mult(size(R, 1), 2, 2) + c_mtx_add(size(y_tilde, 1), size(y_tilde, 2));
    end

    m.s_hat = s_hat;
    m.partb_cmp_arithmetics = partb_cmp_arithmetics;
end