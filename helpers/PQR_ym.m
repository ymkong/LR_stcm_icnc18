%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong18apply] Applying Lattice Reduction Technique to Space-Time Coded Multiplexing Systems 
% This program calls a generalized pairwise QR deomposition, that works together with pairwise ELR, and records the complexity of the procedure
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q, R, p, c] = PQR_ym(H, is_sorted, det_type)

% sorted, ZF
% sorted, MMSE
% un-sorted, ZF
% un-sorted, MMSE

[m, n] = size(H);
p = 1 : 1 : n;

if is_sorted == 1 && strcmp(det_type, 'ZF')
	[Q, R, p] = sorted_pairwiseQR(H, 'H');
	c = m * n^2 - 1/2* m * n; % record complexity

elseif is_sorted == 1 && strcmp(det_type, 'MMSE')
	[Q, R, p] = sorted_pairwiseQR(H, 'Hbar');
	c = m * n^2 - 1/2* m * n; % record complexity

elseif is_sorted == 0 && strcmp(det_type, 'ZF')
	[Q, R] = pairwiseQR(H, 'H');
	c = m * n^2 - 1/2* m * n; % record complexity

elseif is_sorted == 0 && strcmp(det_type, 'MMSE')
	[Q, R] = pairwiseQR(H, 'Hbar');
	c = m * n^2 - 1/2* m * n; % record complexity
end

end