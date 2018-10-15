%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reference paper: [kong18apply] Applying Lattice Reduction Technique to Space-Time Coded Multiplexing Systems 
% This program compares BER and complexity of various direct- and group- detectors in a space-time coded multiplexing systems specified in paper [kong18apply]
% 
% Written by: Yiming Kong
% Date: 3/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc; 
addpath('./det_group', './det_direct', './lr_m', './core_c','./helpers', './mtx_complexity');

% setup simulation model
model = struct(...
    'sim_type', 'SNR_VS_BER', ... % NT_VS_BER, SNR_VS_BER
	'seed', 16, ...
    'STCtype', 'alamouti', ... 
    'sim_n', 1e5, ...
    'chn_n', 1e0, ...
    'algs', [], ...
    'hmod', modem.qammod('M', 16, 'SymbolOrder', 'gray', 'InputType', 'bit'), ...
    'max_no_error', 1e3);

algs = {};
%% GROUP detectors
% % % if strcmp(model.STCtype, 'alamouti')
% % % %     algs{length(algs) + 1} = struct('sn', 'g-zf', 'title', 'Group ZF', 'marker', 'k*-', 'func', @(act, mdl, m) det_group_zf(act, mdl, m, @(H) eyeH(H), 'no_lr_g_zf'));
% % %     algs{length(algs) + 1} = struct('sn', 'g-mmse', 'title', 'Group MMSE', 'marker', 'b>-', 'func', @(act, mdl, m) det_group_mmse(act, mdl, m, @(H) eyeH(H), 'no_lr_g_mmse'));
% % %     algs{length(algs) + 1} = struct('sn', 'g-zf-sic', 'title', 'Group ZF-SIC', 'marker', 'r^-', 'func', @(act, mdl, m) det_group_zf_sic(act, mdl, m, @(H) eyeH(H), 'no_lr_g_zfsic'));
% % % end

%% Direct detectors
%% SIC
% algs{length(algs) + 1} = struct('sn', 'mmse-sic', 'title', 'Direct MMSE-SIC', 'marker', 'm*-', 'func', @(act, mdl, m) det_direct_lr_mmse_sic_generic(act, mdl, m, @(H) eyeH(H), 'mmsesic'));

%% LR-aided SIC
algs{length(algs) + 1} = struct('sn', 'd-clll-mmse-sic', 'title', 'Direct CLLL-MMSE-SIC-0.75', 'marker', 'r*-', 'func', @(act, mdl, m) det_direct_lr_mmse_sic_generic(act, mdl, m, @(H) clll_c(H, 0.75), 'd_clll_mmsesic'));
% algs{length(algs) + 1} = struct('sn', 'd-elr-slb-mmse-sic', 'title', 'Direct ELR-SLB-MMSE-SIC', 'marker', 'm<-', 'func', @(act, mdl, m) det_direct_lr_mmse_sic_generic(act, mdl, m, @(H) elr_dual_c(H, 'elr_slb'), 'elr_slb_dual_mmsesic'));
% algs{length(algs) + 1} = struct('sn', 'd-elr-mmse-sic', 'title', 'Direct ELR-MMSE-SIC', 'marker', 'k<-', 'func', @(act, mdl, m) det_direct_lr_mmse_sic_generic(act, mdl, m, @(H) elr_dual_c(H, 'elr'), 'elr_dual_mmsesic'));
% algs{length(algs) + 1} = struct('sn', 'd-pelrp-mmse-sic', 'title', 'Direct PELR-MMSE-SIC', 'marker', 'g<-', 'func', @(act, mdl, m) det_direct_lr_mmse_sic_generic(act, mdl, m, @(H) elr_dual_c(H, 'pelrp'), 'pelrp_dual_mmsesic'));
% algs{length(algs) + 1} = struct('sn', 'd-pelrp-slb-mmse-sic', 'title', 'Direct PELRP-SLB-MMSE-SIC-YM', 'marker', 'r<-', 'func', @(act, mdl, m) det_direct_lr_mmse_sic_generic(act, mdl, m, @(H) elr_dual_c(H, 'pelrp_slb_ym'), 'pelrp_slb_dual_mmsesic_ym'));

% add detection algorithms
model.algs = algs;

% setup legend for plotting
le = cell(length(algs), 1);
for alg_ind = 1 : length(algs)
    le{alg_ind} = algs{alg_ind}.title; 
end

% based on simulation setup
if strcmp(model.sim_type, 'SNR_VS_BER') % run SNR vs. BER curve, i.e., Fig. 2 in [kong18apply]
    
    disp('running for SNR_VS_BER');

    % add model parameters
    model.Nt = 80;
    model.Nr = 40;
    model.SNRdb = 20 : 5 : 40;

    % run simulator
    [R] = stsm_simulator(model);

    % plot results
    fields = fieldnames(R);
    for k = 1 : numel(fields)
    	if (~strcmp(fields{k}, 'Pe'))
    		figure
		    for alg_ind = 1 : length(algs)
		        plot(model.SNRdb, R.(fields{k})(alg_ind, :), algs{alg_ind}.marker, 'LineWidth', 2); hold on
		    end
		    xlabel('SNR(dB)'); ylabel(fields{k}); legend(le); 
		else
		    figure
		    for alg_ind = 1 : length(algs)
		        semilogy(model.SNRdb, R.(fields{k})(alg_ind, :), algs{alg_ind}.marker, 'LineWidth', 2); hold on
		    end
		    xlabel('SNR(dB)'); ylabel(fields{k}); legend(le);
		end
    end

elseif strcmp(model.sim_type, 'NT_VS_BER') % run Nt vs. SNR curve, i.e., Fig. 1 in [kong18apply]

    disp('running for NT_VS_BER');
    
    % add model parameters
    Nt = 20 : 20 : 100;
    Nr = Nt/2;
    model.SNRdb = 25;
    
    % set up struct variable to collect results
    RR = struct('Pe', [], 'iterations', [], 'basis_updates', [], 'cmp_arithmetics', [], 'arithmetics', [], 'det_cmp_arithmetics', []);

    for k = 1 : length(Nt)
        model.Nr = Nr(k);
        model.Nt = Nt(k);
        
        % run simulator
        [R] = stsm_simulator(model);

        fields = fieldnames(R);
        for f_id = 1 : numel(fields)
            RR.(fields{f_id}) = [RR.(fields{f_id}), R.(fields{f_id})];
        end
    end

    % plot results
    for k = 1 : numel(fields)
        if ~strcmp(fields{k}, 'Pe')
            figure
            for alg_ind = 1 : length(algs)
                plot(Nt, RR.(fields{k})(alg_ind, :), algs{alg_ind}.marker, 'LineWidth', 2); hold on
            end
            xlabel('Nt'); ylabel(fields{k}); legend(le); 
        else
            figure
            for alg_ind = 1 : length(algs)
                semilogy(Nt, RR.(fields{k})(alg_ind, :), algs{alg_ind}.marker, 'LineWidth', 2); hold on
            end
            xlabel('Nt'); ylabel(fields{k}); legend(le);
        end
    end

else 
    disp('simulation type unknown.');
end