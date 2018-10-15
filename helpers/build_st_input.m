%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong18apply] Applying Lattice Reduction Technique to Space-Time Coded Multiplexing Systems 
% This program builds the input signal matrix for the group-detection model
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s_ST = build_st_input(s, STCtype)

Nt = length(s);

if strcmp(STCtype, 'alamouti')

    s2 = [];
    for k = 1 : Nt / 2
        s2 = [s2, s(2*k), -s(2*k-1)];
    end
    s_ST = [s, s2'];
else
    disp('ERROR: unknown type of space-time code');
end

end