%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong18apply] Applying Lattice Reduction Technique to Space-Time Coded Multiplexing Systems 
% This program builds the input signal vector for the direct detection model
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s_d = build_direct_input(s, STCtype)

Nt = length(s);

if strcmp(STCtype, 'alamouti')
    s_d = s([1:2:Nt, 2:2:Nt], :);
else
    disp('ERROR: unknown type of space-time code');
end
    
end