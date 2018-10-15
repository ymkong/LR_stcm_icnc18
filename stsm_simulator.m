%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong18apply] Applying Lattice Reduction Technique to Space-Time Coded Multiplexing Systems 
% This program runs simulation based on parameters specified in the model struct
% 
% Written by: Yiming Kong
% Date: 3/1/2017
% Acknowledgement: the code style in this simulator is heavily influenced by that of Qi Zhou. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results] = stsm_simulator(model)

% set random seed
rand('state', model.seed)
randn('state', model.seed)

% copy down the parameters in the model
sim_n = model.sim_n;
chn_n = model.chn_n;
Nr = model.Nr;
Nt = model.Nt;
SNRdb = model.SNRdb;
algs = model.algs;
hmod = model.hmod;
hdemod = modem.qamdemod(hmod, 'OutputType', 'bit');
Es = mean(abs(hmod.Constellation).^2);
bpers = log2(hmod.M);
Eb = Es / bpers;
sigmas = sqrt(Nt * Eb ./ (10 .^ (SNRdb ./ 10)));
model.sigmas = sigmas;

% setup detector model
det_models = cell(length(algs), 1);
for alg_ind = 1 : length(algs)
    det_models{alg_ind} = struct(...
        'y', [], 'simga', [], 'SNR_ind', 0, ...
        'Es', Es, 'Eb', Eb, 'bpers', bpers, ...
        's', [], 's_d', [], ...
        's_hat', [], ...
        'iterations', [], ...
        'basis_updates', [], ...
        'cmp_arithmetics', [], ...
        'arithmetics', [], ...
        'extra', [], ...
        'partb_cmp_arithmetics', []);
end

% variables to keep track of number of signal blocks run for each SNR
% when the number of errors collected a certain SNR exceeds model.max_no_error, we think that the number of simulations for that SNR is enough
start_SNR_ind = ones(length(algs), 1);
blocks_per_SNR = zeros(length(algs), length(SNRdb));

% variables to record BER and number of arithmetic operations
Pe = zeros(length(algs), length(SNRdb));
iterations = zeros(length(algs), length(SNRdb));
basis_updates = zeros(length(algs), length(SNRdb));
cmp_arithmetics = zeros(length(algs), length(SNRdb));
arithmetics = zeros(length(algs), length(SNRdb));
det_cmp_arithmetics = zeros(length(algs), length(SNRdb)); % number of arithmetric operations for the full detector

tic
for runs = 1 : sim_n % start simulation    
    H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2); % (Nr, Nt)
    
    H_p = build_direct_chn_pair(H, model.STCtype); % equivalent PAIR-WISE direct channel matrix: (2Nr, Nt), or (4Nr, Nt)
    
    model.chn_info_st = struct();
    model.chn_info_di = struct();    
    model.chn_info_st.H = H;
    model.chn_info_di.H_p = H_p; 
    
    for alg_ind = 1 : length(algs)
        
        % run preprocessing steps
        [model, det_models{alg_ind}] = algs{alg_ind}.func('updateH', model, det_models{alg_ind});

        % record the complexity of preprocessing steps
        if (size(det_models{alg_ind}.iterations, 1) * size(det_models{alg_ind}.iterations, 2) > 1) % update results depend on SNR
        	for SNR_ind = 1 : length(SNRdb)
                
                if (start_SNR_ind(alg_ind) > SNR_ind)
                    continue
                end
                
                % update complexity results
        		iterations(alg_ind, SNR_ind) = iterations(alg_ind, SNR_ind) + det_models{alg_ind}.iterations(SNR_ind);
        		basis_updates(alg_ind, SNR_ind) = basis_updates(alg_ind, SNR_ind) + det_models{alg_ind}.basis_updates(SNR_ind);
        		cmp_arithmetics(alg_ind, SNR_ind) = cmp_arithmetics(alg_ind, SNR_ind) + det_models{alg_ind}.cmp_arithmetics(SNR_ind);
        		arithmetics(alg_ind, SNR_ind) = arithmetics(alg_ind, SNR_ind) + det_models{alg_ind}.arithmetics(SNR_ind);
                det_cmp_arithmetics(alg_ind, SNR_ind) = det_cmp_arithmetics(alg_ind, SNR_ind) + det_models{alg_ind}.cmp_arithmetics(SNR_ind) + det_models{alg_ind}.extra(SNR_ind);
        	end
        else
        	for SNR_ind = 1 : length(SNRdb)
                
                if (start_SNR_ind(alg_ind) > SNR_ind)
                    continue
                end
                
                % update complexity results                
        		iterations(alg_ind, SNR_ind) = iterations(alg_ind, SNR_ind) + det_models{alg_ind}.iterations;
        		basis_updates(alg_ind, SNR_ind) = basis_updates(alg_ind, SNR_ind) + det_models{alg_ind}.basis_updates;
        		cmp_arithmetics(alg_ind, SNR_ind) = cmp_arithmetics(alg_ind, SNR_ind) + det_models{alg_ind}.cmp_arithmetics;
        		arithmetics(alg_ind, SNR_ind) = arithmetics(alg_ind, SNR_ind) + det_models{alg_ind}.arithmetics;
                det_cmp_arithmetics(alg_ind, SNR_ind) = det_cmp_arithmetics(alg_ind, SNR_ind) + det_models{alg_ind}.cmp_arithmetics + det_models{alg_ind}.extra;
        	end
        end
    end

    % run simulation
    for chn_ind = 1 : chn_n
        
        % generates random bits
        b = round(rand(Nt * bpers, 1)); % (Nt*bpers, 1)
        
        % modulate to symbols
        s = modulate(hmod, b); % (Nt, 1)

        % build the space-time structure
        s_ST = build_st_input(s, model.STCtype); % (Nt, 2) or (Nt, 4)
        
        % build the equivalent input for a direct detector model
        s_d = build_direct_input(s, model.STCtype); % (Nt, 1)
        
        for alg_ind = 1 : length(algs)
            det_models{alg_ind}.s = s;
            det_models{alg_ind}.s_d = s_d;
        end
        
        % go through channel
        y0 = H * s_ST;
        
        % generate noise
        n0 = 1 / sqrt(2) * (randn(Nr, size(s_ST, 2)) + 1j * randn(Nr, size(s_ST, 2)));
        
        for SNR_ind = 1 : length(SNRdb)
            y = y0 + sigmas(SNR_ind) * n0; % (Nr, 2) or (Nr, 4)

            for alg_ind = 1 : length(algs)

                if (start_SNR_ind(alg_ind) > SNR_ind)
                    continue
                else
                    blocks_per_SNR(alg_ind, SNR_ind) = blocks_per_SNR(alg_ind, SNR_ind) + 1;
                end

                % setup neccesary values for the detector
                det_models{alg_ind}.y = y;
                det_models{alg_ind}.sigma = sigmas(SNR_ind);
                det_models{alg_ind}.SNR_ind = SNR_ind;
                
                % detection
                [model, det_models{alg_ind}] = algs{alg_ind}.func('det', model, det_models{alg_ind});

                % record BER and detection complexity
                Pe(alg_ind, SNR_ind) = Pe(alg_ind, SNR_ind) + sum(demodulate(hdemod, det_models{alg_ind}.s_hat) ~= b);
                det_cmp_arithmetics(alg_ind, SNR_ind) = det_cmp_arithmetics(alg_ind, SNR_ind) + det_models{alg_ind}.partb_cmp_arithmetics;
            end
        end

        % update starting SNR if the number of errors collected for certain SNR is large enough
        for alg_ind = 1 : length(algs)
            while (start_SNR_ind(alg_ind) <= length(SNRdb)) && (Pe(alg_ind, start_SNR_ind(alg_ind)) > model.max_no_error)
                start_SNR_ind(alg_ind) = start_SNR_ind(alg_ind) + 1;
            end
        end
    end

    % print out left running time
    if(mod(runs, 1000) == 0)
        fprintf('iter: %d, left time: %0.2f\n', runs, toc / runs * (sim_n - runs));
    end
end
toc

% compute results
Pe = Pe ./ blocks_per_SNR ./ Nt ./ bpers;
det_cmp_arithmetics = det_cmp_arithmetics ./ blocks_per_SNR ./ Nt ./ bpers;
iterations = iterations ./ (blocks_per_SNR ./chn_n);
basis_updates = basis_updates ./ (blocks_per_SNR ./chn_n);
cmp_arithmetics = cmp_arithmetics ./ (blocks_per_SNR ./chn_n);
arithmetics = arithmetics ./ (blocks_per_SNR ./chn_n);

% tie them to a nob
results = struct('Pe', Pe, ...
	'iterations', iterations, ...
	'basis_updates', basis_updates, ...
	'cmp_arithmetics', cmp_arithmetics, ...
	'arithmetics', arithmetics, ...
	'det_cmp_arithmetics', det_cmp_arithmetics);