%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong18apply] Applying Lattice Reduction Technique to Space-Time Coded Multiplexing Systems 
% This program runs a generic MMSE-SIC based detector given the LR function, and records the complexity of the procedure
% 
% Written by: Yiming Kong
% Date: 3/1/2017
% Acknowledgement: the code style in this simulator is heavily influenced by that of Qi Zhou. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ mdl, m ] = det_direct_lr_mmse_sic_generic(act, mdl, m, lr_func, lr_name)

if (strcmp(act, 'updateH'))

    % get channel matrix
    Hf = mdl.chn_info_di.H_p;
    
    % setup variables
    [mm, nn] = size(Hf);
    NS = length(mdl.SNRdb);
    Ts = zeros(nn, nn, NS);
    ITs = zeros(nn, nn, NS);
    Hs = zeros(mm + nn, nn, NS);
    Qts = zeros(mm + nn, nn, NS);
    Rts = zeros(nn, nn, NS);

    % set up counters to record complexity
    iterations = zeros(1, NS);
    basis_updates = zeros(1, NS);
    cmp_arithmetics = zeros(1, NS);
    arithmetics = zeros(1, NS);
    extra = zeros(1, NS);
    
    for s_ind = 1 : NS
        sigma = mdl.sigmas(s_ind);
        H = [Hf; eye(nn) .* sigma ./ sqrt(m.Es)];
        
        % lattice reduction
        [H_t, T, info] = lr_func(H);

        if (sum(strfind(lr_name, 'pelrp')) > 0) % pairwise elr
            [Qt, Rt, ~, c] = PQR_ym(H_t, 0, 'MMSE');
        elseif (sum(strfind(lr_name, 'clll')) > 0) % clll
            [Qt, Rt, ~, ~] = QR_ym(H_t, 0);
            c = 0;
        else % others
            [Qt, Rt, ~, c] = QR_ym(H_t, 0);
        end

        % record matrices
        Ts(: , : , s_ind) = T;
        ITs(: , : , s_ind) = inv(T);
        Hs(: , : , s_ind) = H_t;
        Qts(: , : , s_ind) = Qt;
        Rts(: , : , s_ind) = Rt;

        % record complexity
        iterations(s_ind) = info.iterations;
        basis_updates(s_ind) = info.basis_updates;
        cmp_arithmetics(s_ind) = info.cmp_arithmetics;
        arithmetics(s_ind) = info.arithmetics;
        extra(s_ind) = c;
    end

    % record channels and supporting matrices
    mdl.chn_info_di.(lr_name).Hs = Hs;
    mdl.chn_info_di.(lr_name).Ts = Ts;
    mdl.chn_info_di.(lr_name).ITs = ITs;
    mdl.chn_info_di.(lr_name).Qts = Qts;
    mdl.chn_info_di.(lr_name).Rts = Rts;
    
    % record complexity
    m.iterations = iterations;
    m.basis_updates = basis_updates;
    m.cmp_arithmetics = cmp_arithmetics;
    m.arithmetics = arithmetics;
    m.extra = extra;

elseif (strcmp(act, 'det'))

    Nt = mdl.Nt;
    partb_cmp_arithmetics = 0;

    y = [];
    if (size(m.y, 2) == 2)  % alamouti code
        y = [m.y(:, 1); conj(m.y(:, 2)); zeros(Nt, 1)];
    end
    
    % get supporting matrices
    Qt = mdl.chn_info_di.(lr_name).Qts(: , : , m.SNR_ind);
    Rt = mdl.chn_info_di.(lr_name).Rts(: , : , m.SNR_ind);
    TC = mdl.chn_info_di.(lr_name).Ts(: , : ,  m.SNR_ind);
    ITC = mdl.chn_info_di.(lr_name).ITs(: , : ,  m.SNR_ind);

    x_nc = Qt' * y/2 - Rt * ITC * 0.5 * (1 + 1j) * ones(Nt, 1);
    partb_cmp_arithmetics = partb_cmp_arithmetics + 2 * size(Qt, 2) * size(Qt, 1); % count complexity

    % SIC detecting
    z_hat = zeros(Nt, 1);
    for nn = Nt : (-1) : 1
        z_hat(nn) = round(x_nc(nn) / Rt(nn, nn));
        x_nc = x_nc - z_hat(nn) * Rt(:, nn);
        partb_cmp_arithmetics = partb_cmp_arithmetics + 1 + 2 * nn; % count complexity
    end

    % transfer to the original domain
    zc_hat = round(TC * z_hat);    
    s_hat = 2 * (zc_hat + (0.5 + 1j * 0.5));
    partb_cmp_arithmetics = partb_cmp_arithmetics + 2 * size(TC, 1) * size(TC, 2); % count complexity
    
    % return results
    m.s_hat = s_hat;
    m.partb_cmp_arithmetics = partb_cmp_arithmetics;

end

end