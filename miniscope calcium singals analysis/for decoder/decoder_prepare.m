%% 在做decode 前获得神经元每个trial的均值 Y 神经元指定类型； X 用于解码的数据
%%
load('E:\data\Ca_analyzing\Calb1_all_ca\final\dffstructcalb1tail.mat');
load('E:\data\Ca_analyzing\Calb1_all_ca\final\lickactivated');
lick_activated=ROIIDResponsive;
load('E:\data\Ca_analyzing\Calb1_all_ca\final\attackactivated');
attack_activated=ROIIDResponsive;
load('E:\data\Ca_analyzing\Calb1_all_ca\final\pinchactivated');
pinch_activated=ROIIDResponsive;
% load('E:\data\tnc miniscope\Tnc-all\final\dffstructtnctail.mat');
% load('E:\data\tnc miniscope\Tnc-all\final\lickactivated');
% lick_activated=ROIIDResponsive;
% load('E:\data\tnc miniscope\Tnc-all\final\attackactivated');
% attack_activated=ROIIDResponsive;
% load('E:\data\tnc miniscope\Tnc-all\final\pinchactivated');
% pinch_activated=ROIIDResponsive;


activated_neuron=[lick_activated,attack_activated,pinch_activated];
sample_no=length (unique (fieldnames (dffstruct)));
get_dff=dff_data;
[roi_no]=get_dff.getroi_no(dffstruct);
%根据行为筛选
[tail_licking,tail_l_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'l','2');
[tail_attacking,tail_a_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'u','2');%,'2');
[tail_clip,tail_j_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'j','2');
% zscore calculate
Fs=14.854;
baseline=int32(5*Fs);
% cell_means = cellfun(@(x) mean(x(:, 1:baseline), 2), tail_licking, 'UniformOutput', false);
% cell_stds = cellfun(@(x) std(x(:, 1:baseline), 0, 2), tail_licking, 'UniformOutput', false);
% lick_ca_zscore = cellfun(@(x, m, s) (x - m) ./ s,tail_licking, cell_means, cell_stds, 'UniformOutput', false);
% 
% cell_means = cellfun(@(x) mean(x(:, 1:baseline), 2),tail_attacking, 'UniformOutput', false);
% cell_stds = cellfun(@(x) std(x(:, 1:baseline), 0, 2),tail_attacking, 'UniformOutput', false);
% attack_ca_zscore = cellfun(@(x, m, s) (x - m) ./ s,tail_attacking, cell_means, cell_stds, 'UniformOutput', false);
% 
% cell_means = cellfun(@(x) mean(x(:, 1:baseline), 2),tail_clip, 'UniformOutput', false);
% cell_stds = cellfun(@(x) std(x(:, 1:baseline), 0, 2),tail_clip, 'UniformOutput', false);
% pinch_ca_zscore = cellfun(@(x, m, s) (x - m) ./ s,tail_clip, cell_means, cell_stds, 'UniformOutput', false);

% tail_licking average
fra_leng=size(tail_licking{1},2);
dff_mean = cell(1,roi_no);
for p = 1:roi_no
    dff_mean{p} = mean(tail_licking{p});
end % trial 的平均
lick_av_dff = zeros(roi_no,fra_leng);
for q = 1:roi_no
    lick_av_dff(q,:) = dff_mean{q};
end
% lick_av_dff=lick_av_dff(:,baseline:222);
lick_av_dff = normalize(lick_av_dff, 2, 'range', [0 1]); % 按行归一化到 [0,1]
lick_av_dff = lick_av_dff * 2 - 1; % 变换到 [-1,1]
% save('E:\data\Ca_analyzing\Calb1_all_ca\final\lick_av_dff','lick_av_dff');


% % tail_attack average
fra_leng=size(tail_attacking{1},2);
dff_mean = cell(1,roi_no);
for p = 1:roi_no
    dff_mean{p} = mean(tail_attacking{p});
end % trial 的平均
attack_av_dff = zeros(roi_no,fra_leng);
for q = 1:roi_no
    attack_av_dff(q,:) = dff_mean{q};
end
% attack_av_dff=attack_av_dff(:,baseline:222);
attack_av_dff = normalize(attack_av_dff, 2, 'range', [0 1]); % 按行归一化到 [0,1]
attack_av_dff = attack_av_dff * 2 - 1; % 变换到 [-1,1]
% save('E:\data\Ca_analyzing\Calb1_all_ca\final\attack_av_dff','attack_av_dff');


% % % % tail_clip average
% % % fra_leng=size(tail_clip{1},2);
% % % dff_mean = cell(1,roi_no);
% % % dff_sem = cell(1,roi_no);
% % % for p = 1:roi_no
% % %     dff_mean{p} = mean(tail_clip{p});
% % %     dff_sem{p} = std(tail_clip{p})/sqrt(size(tail_clip{p},1));
% % % end % trial 的平均
% % % pinch_av_dff = zeros(roi_no,fra_leng);
% % % for q = 1:roi_no
% % %     pinch_av_dff(q,:) = dff_mean{q};
% % % end
% % % pinch_av_dff=pinch_av_dff(:,baseline:222);
% % % % save('E:\data\Ca_analyzing\Calb1_all_ca\for PCA\pinch_av_dff','pinch_av_dff');



% 准备 X 文件 for lick veus attack
lick_veus_attack=[lick_av_dff;attack_av_dff];% 固定lick 然后是attack 然后是pinch 注意不要改变，改变后续程序需修改
% zscore_average=[]

%准备Y文件
y_decoder_lick=zeros(1,roi_no);
for roi=1:roi_no
    if ismember(roi,lick_activated)
        y_decoder_lick(roi)=1;%lick 编码
%     elseif ismember(roi,attack_activated)
%         y_decoder(roi)=2;%attack 编码
%     elseif ismember(roi,pinch_activated)
%         y_decoder(roi)=3;%pinch 编码
    end            
end
y_decoder_attack=zeros(1,roi_no);
for roi=1:roi_no
    if ismember(roi,attack_activated)
        y_decoder_attack(roi)=1;%lick 编码
%     elseif ismember(roi,attack_activated)
%         y_decoder(roi)=2;%attack 编码
%     elseif ismember(roi,pinch_activated)
%         y_decoder(roi)=3;%pinch 编码
    end            
end


y_decoder_lickattack=zeros(1,roi_no);
for roi=1:roi_no
    if ismember(roi,lick_activated)
        y_decoder_lickattack(roi)=1;%lick 编码
    elseif ismember(roi,attack_activated)
        y_decoder_lickattack(roi)=2;%attack 编码
%     elseif ismember(roi,pinch_activated)
%         y_decoder(roi)=3;%pinch 编码
    end            
end
% save('E:\data\Ca_analyzing\Calb1_all_ca\for PCA/y_decoder.mat','y_decoder');
% save('E:\data\tnc miniscope\Tnc-all\final\decoder/y_decoder_lick.mat','y_decoder_lick');
save('E:\data\Ca_analyzing\Calb1_all_ca\final\decoder/y_decoder_lickattack.mat','y_decoder_lickattack');
% save('E:\data\tnc miniscope\Tnc-all\final\decoder/y_decoder_attack.mat','y_decoder_attack');
% save('E:\data\Ca_analyzing\Calb1_all_ca\for PCA/dff_average.mat','dff_average');
% save('E:\data\Ca_analyzing\Calb1_all_ca\final\decoder/lick_veus_attack.mat','lick_veus_attack');
%%
tail_lick_trials = [];
for i = 1:length(tail_licking)
    data = tail_licking{i};  % n×222
    tail_lick_trials = [tail_lick_trials; data];
end
%% shuffle data
% 
[lickingevent]=reorganize(tail_licking);
[attackevent]=reorganize(tail_attacking);
[pinchevent]=reorganize(tail_clip);
cellid=1;
lick_shuffle=cell(1,roi_no);
attack_shuffle=cell(1,roi_no);
for i=1:length(lickingevent)
    [trials1,time,neurons]=size(lickingevent{i});
    [trials2,time,neurons]=size(attackevent{i});
    for j=1:neurons
        D_permuation = [lickingevent{i}(:,:,j);attackevent{i}(:,:,j)];
        sel_p=randperm(trials1+trials2);
        D_permuation=D_permuation(sel_p,:);
        D1_perm = D_permuation(1:trials1,:);
        D2_perm = D_permuation(trials1+1:end,:);% 原理就是把event标签打乱
        lick_shuffle{cellid}=D1_perm;
        attack_shuffle{cellid}=D2_perm;
        cellid=cellid+1;
    end
end

fra_leng=size(lick_shuffle{1},2);
dff_mean = cell(1,roi_no);
for p = 1:roi_no
    dff_mean{p} = mean(lick_shuffle{p});
end % trial 的平均
lick_shuffle_av_dff = zeros(roi_no,fra_leng);
for q = 1:roi_no
    lick_shuffle_av_dff(q,:) = dff_mean{q};
end
% lick_shuffle_av_dff=lick_shuffle_av_dff(:,baseline:222);
lick_shuffle_av_dff = normalize(lick_shuffle_av_dff, 2, 'range', [0 1]); % 按行归一化到 [0,1]
lick_shuffle_av_dff = lick_shuffle_av_dff * 2 - 1; % 变换到 [-1,1]
% save('E:\data\Ca_analyzing\Calb1_all_ca\final\decoder/lick_shuffle_av_dff.mat','lick_shuffle_av_dff');

fra_leng=size(attack_shuffle{1},2);
dff_mean = cell(1,roi_no);
for p = 1:roi_no
    dff_mean{p} = mean(attack_shuffle{p});
end % trial 的平均
attack_shuffle_av_dff = zeros(roi_no,fra_leng);
for q = 1:roi_no
    attack_shuffle_av_dff(q,:) = dff_mean{q};
end
% lick_shuffle_av_dff=lick_shuffle_av_dff(:,baseline:222);
attack_shuffle_av_dff = normalize(attack_shuffle_av_dff, 2, 'range', [0 1]); % 按行归一化到 [0,1]
attack_shuffle_av_dff= attack_shuffle_av_dff  * 2 - 1; % 变换到 [-1,1]
% save('E:\data\Ca_analyzing\Calb1_all_ca\final\decoder/attack_shuffle_av_dff.mat','attack_shuffle_av_dff');
% 
twoeventshuffle=[lick_shuffle_av_dff;attack_shuffle_av_dff];
save('E:\data\Ca_analyzing\Calb1_all_ca\final\decoder/twoeventshuffle3.mat','twoeventshuffle')
