%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program builds an identity matrix based on the size of H, and records the complexity of the procedure
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H, T, info] = eyeH(H)

T = eye(size(H , 2));
info = struct('iterations', 0, 'basis_updates', 0, 'cmp_arithmetics', 0, 'arithmetics', 0);

end