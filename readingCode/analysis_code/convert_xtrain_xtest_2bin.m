% This script was used to convert the final states (based on mean firing rates) 
% stored in x_train and x_test (generated by the function make_xtrain_xtest)
% into binaraized versions using a threshold of 30Hz. That is, x_train/x_test 
% contain the mean final firing rates of each unit while x_train_bin/x_test_bin 
% (produced by this script) binarize these responses such that if the mean 
% final firing rate of a unit is >30Hz it has a value of 1 ('on') and if 
% the mean final firing rate is <30Hz it has a value of 0 ('off').
datadirs{1} = 'We0_Wee_sweep_noise000';
datadirs{2} = 'We0_Wee_sweep_noise002';
datadirs{3} = 'dif_start_We0_Wee_sweep_noise000';
datadirs{4} = 'dif_start59_We0_Wee_noise000';
datadirs{5} = 'dif_start59_stim20_We0_Wee_noise000';
datadirs{6} = 'We0_Wee_sweep_noise000_amp100_110';
datadirs{7} = 'no_D_We0_Wee_sweep_noise000';
datadirs{8} = 'Wei_Wie_sweep_noise000';
datadirs{9} = 'Wei_Wie_sweep_noise002';
datadirs{10} = 'variable_amplitudes_Wei_Wie_sweep_noise000';
datadirs{11} = 'variable_durations_Wei_Wie_sweep_noise000';
datadirs{12} = 'We0_sweep_noise000';
datadirs{13} = 'no_EE_We0_sweep_noise000';
datadirs{14} = 'dif_start20_stim20_We0_Wee_noise000';
datadirs{15} = 'We0_Wee_sweep_noise000_amp100_110';
datadirs{16} = 'We0_Wee_dif_stim_int';
datadirs{17} = 'We0_Wee_word_seq';
datadirs{18} = 'no_D_We0_Wee_sweep_dif_area';
datadirs{19} = 'We0_Wee_Iapp_generalization';
datadirs{20} = 'We0_Wee_word_seq_noise002';
datadirs{21} = 'no_D_We0_Wee_word_seq_dif_area';
dirs2analyze = [19 20 21];
thresh = 30;
for i=1:length(datadirs)
    if (~ismember(i,dirs2analyze))
        continue;
    end
    disp(datadirs{i})
    rl = load(['data/' datadirs{i} '/run_log.mat']); rl=rl.run_log;
    if (isfield(rl,'param1_range'))
        nv1 = length(rl.param1_range);
        nv2 = length(rl.param2_range);
        nnets = rl.nnets;
        for j=1:nv1
            for k=1:nv2 
                for l=1:nnets
                    curdir = ['data/' datadirs{i} '/' num2str(j) '/' num2str(k) '/net' num2str(l)];
                    xtr = load([curdir '/x_train.mat']); xtr=xtr.x;
                    xtst = load([curdir '/x_test.mat']); xtst=xtst.x;
                    xtr_bin = (xtr > thresh);
                    xtst_bin = (xtr > thresh);
                    x = xtr_bin;
                    save([curdir '/x_train_bin.mat'],'x','-mat')
                    x = xtst_bin;
                    save([curdir '/x_test_bin.mat'],'x','-mat')
                end
            end
        end
    else
        nv = length(rl.param_range);
        nnets = rl.nnets;
        for j=1:nv
            for k=1:nnets
                curdir = ['data/' datadirs{i} '/' num2str(j) '/net' num2str(k)];
                xtr = load([curdir '/x_train.mat']); xtr=xtr.x;
                xtst = load([curdir '/x_test.mat']); xtst=xtst.x;
                xtr_bin = (xtr > thresh);
                xtst_bin = (xtr > thresh);
                x = xtr_bin;
                save([curdir '/x_train_bin.mat'],'x','-mat')
                x = xtst_bin;
                save([curdir '/x_test_bin.mat'],'x','-mat')
            end
        end
    end
    disp(['done with ' datadirs{i}])
end