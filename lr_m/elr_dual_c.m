%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calls the (pairwise) element-based lattice reduction algorithms, and records the complexity of the procedure
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H_t, T, info] = elr_dual_c(H, strat)

[m, n] = size(H);

G = (H' * H);
C = inv(G);

if(strcmp(strat, 'elr'))
    [T, info] = elr_dual_core_c(C, eye(n), m);

elseif(strcmp(strat, 'elr_slb'))
    [T, info] = elr_slb_dual_core_c(C, eye(n), m);

elseif(strcmp(strat, 'pelrp'))
    [T, info] = pelrp_dual_core_c(C, eye(n), m);

elseif(strcmp(strat, 'pelrp_slb_ym'))
    [T, info] = pelrp_slb_dual_core_c_ym(C, eye(n), m);
end

info.pre_C = C;
info.pre_G = G;
H_t = H * T;

% 1) add pre-processing complexity: (H'H); (H'H)-1 
% 2) [H_tilde, T] are returned by LR so the complexity of computing H_tilde is counted inside LR.
if (strcmp(strat, 'elr') || strcmp(strat, 'elr_slb')) % no pairwise

    info.cmp_arithmetics = info.cmp_arithmetics + (4/3 * n^3 + 2 * n^2 * m);

elseif (sum(strfind(strat, 'pelrb')) > 0 || sum(strfind(strat, 'pelrp')) > 0)    
    
    info.cmp_arithmetics = info.cmp_arithmetics + (2/3 * n^3 + n^2 * m);
    
end

end