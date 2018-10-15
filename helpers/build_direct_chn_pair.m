%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong18apply] Applying Lattice Reduction Technique to Space-Time Coded Multiplexing Systems 
% This program builds the direct channel based on the direct detection model
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H_p = build_direct_chn_pair(H, STCtype)

[~, Nt] = size(H);

if strcmp(STCtype, 'alamouti')

    % example
    % H = [h1, h2, h3, h4]
    % H_d = [ h1,  h2,   h3,   h4 ]
    %       [-h2*, h1*, -h4*,  h3*]

    H_b = [];
    for k = 1 : Nt / 2
        H_b = [H_b, -conj(H(:, 2 * k)), conj(H(:, 2 * k - 1))];
    end

    H_p = [H; H_b];
    
else
    disp('ERROR: unknown type of space-time code');
end


end