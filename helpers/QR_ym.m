%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calls a QR deomposition, and records the complexity of the procedure
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q, R, p, c] = QR_ym(H, isSorted)

[m, n] = size(H);

p = 1 : 1 : n;

if(isSorted)
    [Q, R, p] = sorted_QR(H);
    c = 2 * m * n^2 + m * n; % record complexity
else
    [Q, R] = qr(H, 0);
    c = 2 * m * n^2 + m * n; % record complexity
end