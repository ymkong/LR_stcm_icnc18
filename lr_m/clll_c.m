%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calls the c implementation of CLLL, and records the complexity of the procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ H_t, T, info ] = clll_c( H, info )

if (nargin <= 1)
    info = struct('delta', 0.75);
end

[~, R] = qr(H, 0);
T = eye(size(R, 2));
[T, info] = clll_core_c(R, T, info);
H_t = H * T;

if (nargout >= 2)
    info.R = R;
end

% preprocessing: QR
[m, n] = size(H);
info.cmp_arithmetics = info.cmp_arithmetics + (2 * m * n^2 + m * n);

end