load('E:\data\Ca_analyzing\Calb1_all_ca\final\dffstructcalb1tail.mat')
sample_no=length (unique (fieldnames (dffstruct)));
get_dff=dff_data;
[roi_no]=get_dff.getroi_no(dffstruct);
responseneuron=zeros(sample_no,5);
for i=1:sample_no
    responseneuron(i,1)=i;  
end
%根据行为筛选
% [formalin_dff,formalin_l_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'l','5','2');
% [formalin_flinching,formalin_f_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'f','5');
[tail_licking,tail_l_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'l','2');
% [tail_attacking,tail_a_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'u','2')%,'2');
% [tail_clip,tail_j_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'j','2');
% [tail_control,tail_f_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'f','2');



%%
dff_analyzing=tail_licking; %挑选分析的数据
rawID=tail_l_rawID;
%% 进行基线校正
pre_onset=5;%反应起点(formalin=4, tailpinch=5)
Fs = 14.854;
baseline_range = (floor((pre_onset-5)*Fs)+1):(floor(pre_onset*Fs));
roi_no=size(dff_analyzing,2);
dff_analyzing_1p=cell(1,roi_no);
for rn=1:roi_no
    trace_analyzing=dff_analyzing{rn};
    [n_trials, n_frames] = size(trace_analyzing);
    dff_corrected = zeros(n_trials, n_frames);
    for trial = 1:n_trials
        % 提取当前 trial 的基线范围 ΔF/F 数据
        baseline_dff = trace_analyzing(trial, baseline_range);
        time_baseline = (1:length(baseline_range)) / Fs; % 基线时间点（单位：秒）

        % 线性回归拟合（基线范围）
        p = polyfit(time_baseline, baseline_dff, 1); % 一阶多项式拟合

        % 生成整个 trial 的时间点
        time_full = (1:n_frames) / Fs; 

        % 计算漂移趋势
        trend = polyval(p, time_full); % 漂移趋势（线性拟合值）

        % 扣除漂移趋势，得到校正后的 ΔF/F
        dff_corrected(trial, :) = trace_analyzing(trial, :) - trend;
    end
    dff_analyzing_1p{rn}=dff_corrected;
end

dff_analyzing=dff_analyzing_1p;

%%
fra_leng=size(dff_analyzing{1},2);
dff_mean = cell(1,roi_no);
dff_sem = cell(1,roi_no);
for p = 1:roi_no
    dff_mean{p} = mean(dff_analyzing{p});
    dff_sem{p} = std(dff_analyzing{p})/sqrt(size(dff_analyzing{p},1));
end % trial 的平均
av_dff = zeros(roi_no,fra_leng);
for q = 1:roi_no
    av_dff(q,:) = dff_mean{q}(1:fra_leng);
end

%%  inhibited neurons
pre_onset=5;%反应起点(formalin=4, tailpinch=5)
Fs = 14.854; 
base=2; %挑选baseline 的范围(formalin=0)

ROINuminhibition = 0;
ROIIDinhibiton=[];
average_dff=av_dff;
for rn = 1:roi_no
    TrialNuminhibiton{rn} = 0;
    train_of_appoint= size(dff_analyzing{rn},1);
    response_range = (floor(pre_onset*Fs)+1):(floor((pre_onset+2)*Fs)); % 响应窗口
    base_range=floor(base*Fs)+1:(floor((base+2)*Fs));% baseline 窗口
    for i = 1:train_of_appoint
        if  mean(dff_analyzing{rn}(i,response_range))< mean(dff_analyzing{rn}(i,floor(base*Fs)+1:(floor((base+2)*Fs)))) -1 * std(dff_analyzing{rn}(i,floor(base*Fs)+1:(floor((base+2)*Fs))))   
            TrialNuminhibiton{rn} = TrialNuminhibiton{rn} + 1;
        else
            TrialNuminhibiton{rn} = TrialNuminhibiton{rn};
        end
    end
    ResponseWindow{rn} = average_dff(rn,(floor((pre_onset)*Fs)+1):floor((pre_onset+2)*Fs));
    PostWindow{rn} = average_dff(rn,:);
    min_dff = min(average_dff(rn,response_range));

    if TrialNuminhibiton{rn} > 0.5 * train_of_appoint && mean(ResponseWindow{rn}) <mean(average_dff(rn,base_range))- 2 * std(average_dff(rn,base_range))&& min_dff < mean(average_dff(rn,base_range)) - 3 * std(average_dff(rn,base_range))
            ROINuminhibition = ROINuminhibition + 1;
            ROIIDinhibiton(ROINuminhibition)=rn;
            TrialPerResponsive(ROINuminhibition)=TrialNuminhibiton{rn}/train_of_appoint*100;            
    end
end
%% activated neurons
pre_onset=5;%反应起点(formalin=4, tailpinch=5)
Fs = 14.854; 
base=2; %挑选baseline 的范围(formalin=0)
StimDuration=5;
ROINumresponsive = 0;
ROIIDResponsive=[];
sustained_neun=[];
LatencyToOnset=[];
average_dff=av_dff;
sustained_count=0;
for rn = 1:roi_no
    TrialNumactivation{rn} = 0;
    train_of_appoint= size(dff_analyzing{rn},1);
    response_range = (floor(pre_onset*Fs)+1):(floor((pre_onset+3)*Fs)); % 响应窗口
    base_range=floor(base*Fs)+1:(floor((base+3)*Fs));% baseline 窗口
    for i = 1:train_of_appoint
        if  mean(dff_analyzing{rn}(i,response_range))> 1 * std(dff_analyzing{rn}(i,floor(base*Fs)+1:(floor((base+3)*Fs))))   
            TrialNumactivation{rn} = TrialNumactivation{rn} + 1;
        else
            TrialNumactivation{rn} = TrialNumactivation{rn};
        end
    end
    ResponseWindow{rn} = average_dff(rn,(floor((pre_onset)*Fs)+1):floor((pre_onset+3)*Fs));
    PostWindow{rn} = average_dff(rn,:);
    max_dff = max(average_dff(rn,response_range));

    if TrialNumactivation{rn} >= 0.5 * train_of_appoint && mean(ResponseWindow{rn}) > 2 * std(average_dff(rn,base_range))&& max_dff >  3* std(average_dff(rn,base_range))
            ROINumresponsive = ROINumresponsive + 1;            
            ROIIDResponsive(ROINumresponsive)=rn;
            TrialPerResponsive(ROINumresponsive)=TrialNumactivation{rn}/train_of_appoint*100;
            AmplitudeOfPeak(ROINumresponsive)=max(PostWindow{rn});
            PeakPoint(ROINumresponsive)=find(PostWindow{rn}==AmplitudeOfPeak(ROINumresponsive));
            OnsetThreshold(ROINumresponsive)=find(PostWindow{rn}(1:PeakPoint(ROINumresponsive)) > 2 * std(average_dff(rn,(floor(base*Fs)+1):floor((base+3)*Fs))),1);
            for i = OnsetThreshold(ROINumresponsive):PeakPoint(ROINumresponsive)
                OnsetCriteria = find(PostWindow{rn}(i:i+10) <= 2 * std(average_dff(rn,(floor(base*Fs)+1):floor((base+3)*Fs))));
                if isempty(OnsetCriteria)
                    OnsetPoint(ROINumresponsive) = i;
                    LatencyToOnset(ROINumresponsive)=OnsetPoint(ROINumresponsive)/Fs - 5;
                    break;
                end   
                LatencyToOnset(ROINumresponsive) = 9999;
            end 
        
            stim_end_idx = PeakPoint(ROINumresponsive) + round(StimDuration * Fs);  % 从 peak 往后 StimDuration 秒
            stim_end_idx = min(stim_end_idx, length(average_dff(rn,:)));
            post_peak_signal = average_dff(rn,PeakPoint(ROINumresponsive):stim_end_idx);
            mean_after_peak = mean(post_peak_signal);
            ratio(ROINumresponsive) = mean_after_peak / AmplitudeOfPeak(ROINumresponsive); %大于75 为sustained
            if ratio(ROINumresponsive) >= 0.75
               sustained_count = sustained_count + 1;
               sustained_neun=[sustained_neun,rn];
            end
    end
end
figure;
edges = -1:0.2:2;
histogram(LatencyToOnset, edges, 'Normalization', 'probability', 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
hold on
plot([0,0],[0,0.3],'b--','LineWidth',2);
yticks(0:0.1:0.3);  

if ratio >= 0.75
        sustained_count = sustained_count + 1;
end

%%
dff_analyzing=tail_clip; %挑选分析的数据
rawID=tail_j_rawID;
fra_leng=size(dff_analyzing{1},2);
dff_mean = cell(1,roi_no);
dff_sem = cell(1,roi_no);
for p = 1:roi_no
    dff_mean{p} = mean(dff_analyzing{p});
    dff_sem{p} = std(dff_analyzing{p})/sqrt(size(dff_analyzing{p},1));
end % trial 的平均
av_dff = zeros(roi_no,fra_leng);
for q = 1:roi_no
    av_dff(q,:) = dff_mean{q}(1:fra_leng);
end

%%
neuronID=[];
neuronID(:,1)=str2num(rawID(:,4:5));
neuronID(:,2)=str2num(rawID(:,8:9));

%%
onset=5;
dff_response = cell(1,length(ROIIDResponsive));
dff_inhibtion = cell(1,length(ROIIDinhibiton));
dff_noresponse=cell(1,roi_no-length(ROIIDResponsive)-length(ROIIDinhibiton));


r=1;
nr=1;
ir=1;
for i=1:roi_no
    if ismember(i, ROIIDResponsive)==1
        dff_response{r}=dff_analyzing{i};
        r=r+1;
    elseif ismember(i,ROIIDinhibiton)==1
        dff_inhibtion{ir}=dff_analyzing{i};
        ir=ir+1;
    else
        dff_noresponse{nr}=dff_analyzing{i};
        nr=nr+1;
    end
end

dff_response_ppp=[];
for i=1:size(dff_response,2) 
    dff=dff_response{i};
    dff_response_ppp=cat(1,dff_response_ppp,dff);       
end
dff_response_mean=mean(dff_response_ppp);
dff_response_sem=std(dff_response_ppp)/sqrt(size(dff_response_ppp,1));

noresponse=[];
for i=1:size(dff_noresponse,2) 
    dff=dff_noresponse{i};
    noresponse=cat(1,noresponse,dff);       
end
noresponse_mean=mean(noresponse);
noresponse_sem=std(noresponse)/sqrt(size(noresponse,1));

inhibit_trace=[];
for i=1:size(dff_inhibtion,2) 
    dff=dff_inhibtion{i};
    inhibit_trace=cat(1,inhibit_trace,dff);       
end
inhibit_mean=mean(inhibit_trace);
inhibit_sem=std(inhibit_trace)/sqrt(size(inhibit_trace,1));

figure;
subplot(2,4,2)
yMax=max(dff_response_mean + dff_response_sem); yMin=min(dff_response_mean - dff_response_sem);
imagesc(dff_response_ppp);
hold on
plot([onset*Fs+1,onset*Fs+1],[0,size(dff_response_ppp,1)+0.5],'w--','LineWidth',2);
set(gca,'XTickLabel',{'-5','0','10'},...
    'XTick',[1 onset*Fs+1 size(dff_response_ppp,2)]);
set(gca,'tickdir','out');  
set(gca,'Ytick',[],'YtickLabel',[]);
xlabel('Time from action onset (s)','fontsize',30);
ylabel('Neurons x Trials','fontsize',30);
set(gca,'Fontsize',30,'YDir','normal');
box off
caxis([yMin,yMax]);
actucolorbar = colorbar;
set(actucolorbar,'position',[0.51 0.585 0.006 0.15]);
title(actucolorbar,['         \Delta' 'F/F0(%)'],'fontsize',30);
set(gca,'LineWidth',2);
set(gca,'FontSize',30);


subplot(2,4,4)
x=1:size(dff_response_ppp,2);
shadedErrorBar(x,dff_response_mean,dff_response_sem,{'-r','linewidth',2});
hold on
shadedErrorBar(x,noresponse_mean,noresponse_sem,{'-k','linewidth',2});
hold on


plot([onset*Fs+1,onset*Fs+1],[yMin-0.2*abs(yMin),yMax+0.2*abs(yMax)],'k--','LineWidth',2);
set(gca,'XTickLabel',{'-5','0','10'},'XTick',[1 onset*Fs+1 size(dff_response_ppp,2)]);
set(gca,'tickdir','out');     
xlim([0,size(dff_response_ppp,2)]);
ylim([yMin-0.2*abs(yMin),yMax+0.2*abs(yMax)]);
xlabel('Time from pinch onset (s)','fontsize',30);
ylabel(['\Delta' 'F/F0 (%)'],'fontsize',30);
box off
set(gca,'FontSize',30);
set(gca,'LineWidth',2);
%%
subplot(2,4,6)
yMax=max(inhibit_mean + inhibit_sem); yMin=min(inhibit_mean - inhibit_sem);
hold on
plot([onset*Fs+1,onset*Fs+1],[0,size(dff_response_ppp,1)+0.5],'w--','LineWidth',2);
set(gca,'XTickLabel',{'-5','0','10'},...
    'XTick',[1 onset*Fs+1 size(dff_response_ppp,2)]);
set(gca,'tickdir','out');  
set(gca,'Ytick',[],'YtickLabel',[]);
xlabel('Time from action onset (s)','fontsize',30);
ylabel('Neurons x Trials','fontsize',30);
set(gca,'Fontsize',30,'YDir','normal');
box off
caxis([yMin,yMax]);
actucolorbar = colorbar;
set(actucolorbar,'position',[0.51 0.2 0.006 0.15]);
title(actucolorbar,['         \Delta' 'F/F0(%)'],'fontsize',30);
set(gca,'LineWidth',2);
set(gca,'FontSize',30);


subplot(2,4,8)
x=1:size(dff_response_ppp,2);
shadedErrorBar(x,noresponse_mean,noresponse_sem,{'-k','linewidth',2});
hold on
shadedErrorBar(x,inhibit_mean,inhibit_sem,{'-b','linewidth',2});
hold on

plot([onset*Fs+1,onset*Fs+1],[yMin-0.2*abs(yMin),yMax+0.2*abs(yMax)],'k--','LineWidth',2);
set(gca,'XTickLabel',{'-5','0','10'},'XTick',[1 onset*Fs+1 size(dff_response_ppp,2)]);
set(gca,'tickdir','out');     
xlim([0,size(dff_response_ppp,2)]);
ylim([yMin-0.2*abs(yMin),yMax+0.2*abs(yMax)]);
xlabel('Time from pinch onset (s)','fontsize',30);
box off
set(gca,'FontSize',30);
set(gca,'LineWidth',2);





%%
dff_plot=[];
for i=1:size(dff_analyzing,2) 
    dff=mean(dff_analyzing{i},1);
    dff_plot=cat(1,dff_plot,dff);       
end
dff_plot_mean=mean(dff_plot);
dff_plot_sem=std(dff_plot)/sqrt(size(dff_plot,1));

figure;
subplot(2,4,2)
yMax=max(dff_plot_mean + dff_plot_sem); yMin=min(dff_plot_mean - dff_plot_sem);
M=max(dff_plot(:)); N=min(dff_plot(:));
[dff_sorted_,Index] = sortrows(mean(dff_plot(:,(onset-1)*Fs:(onset+2)*Fs),2),-1);
dff_plot_sorted = dff_plot(Index,:);

imagesc(dff_plot_sorted, [N,M]);
hold on
plot([onset*Fs+1,onset*Fs+1],[0,size(dff_plot,1)+0.5],'w--','LineWidth',2);
set(gca,'XTickLabel',{'-5','0','10'},...
    'XTick',[1 onset*Fs+1 size(dff_plot,2)]);
set(gca,'tickdir','out');  
% set(gca,'Ytick',[],'YtickLabel',[]);
xlabel('Time from pinch onset (s)','fontsize',30);
ylabel('Neurons','fontsize',30);
set(gca,'Fontsize',30,'YDir','normal');
box off
caxis([0,20]);
actucolorbar = colorbar;
set(actucolorbar,'position',[0.51 0.7 0.006 0.15]);
title(actucolorbar,['         \Delta' 'F/F0(%)'],'fontsize',25);
set(gca,'LineWidth',2);
set(gca,'FontSize',30);
%saveas(gcf,['D:\图片\23\23.12\Ca\' 'formalin_neurons'  '.png']);


%% 获取神经元钙信号的zscore
baseline=int32(5*Fs);
cell_means = cellfun(@(x) mean(x(:, 1:baseline), 2), dff_analyzing, 'UniformOutput', false);
cell_stds = cellfun(@(x) std(x(:, 1:baseline), 0, 2), dff_analyzing, 'UniformOutput', false);

ca_zscore = cellfun(@(x, m, s) (x - m) ./ s, dff_analyzing, cell_means, cell_stds, 'UniformOutput', false);


%% zscore 热图
zscore_plot=[];
for i=1:size(ca_zscore,2) 
    dff=mean(ca_zscore{i},1);
    zscore_plot=cat(1,zscore_plot,dff);       
end
zscore_plot_mean=mean(zscore_plot);
zscore_plot_sem=std(zscore_plot)/sqrt(size(zscore_plot,1));
% 创建一个新的figure窗口
figure;
subplot(2,4,2)
yMax=max(zscore_plot_mean + zscore_plot_sem); yMin=min(zscore_plot_mean - zscore_plot_sem);
M=max(zscore_plot(:)); N=min(zscore_plot(:));
[zscore_sorted_,Index] = sortrows(mean(zscore_plot(:,(onset-1)*Fs:(onset+2)*Fs),2),-1);
zscore_plot_sorted = zscore_plot(Index,:);
% 使用imagesc函数绘制热图
imagesc(zscore_plot_sorted, [N,M]);
hold on
plot([onset*Fs+1,onset*Fs+1],[0,size(zscore_plot,1)+0.5],'w--','LineWidth',2);
set(gca,'XTickLabel',{'-5','0','10'},...
    'XTick',[1 onset*Fs+1 size(zscore_plot,2)]);
set(gca,'tickdir','out');  
% set(gca,'Ytick',[],'YtickLabel',[]);
xlabel('Time from action onset (s)','fontsize',30);
ylabel('Neurons','fontsize',30);
set(gca,'Fontsize',30,'YDir','normal');
box off


% 定义颜色点
cmap_points = [0.41, 0.35, 0.8;    % 蓝色
               1 0.9 0.85;    % 白色
               0.86 0.07 0.23]; % 粉红色

% 生成平滑的颜色映射
n = 256; % 颜色映射的分辨率
cmap = interp1([0, 0.5, 1], cmap_points, linspace(0, 1, n));

% 应用自定义颜色映射
colormap(cmap);
% 调整颜色范围
caxis([-3, 3]);
actucolorbar = colorbar;
set(actucolorbar,'position',[0.51 0.7 0.006 0.15]);
title(actucolorbar,['  zscore'],'fontsize',25);
set(gca,'LineWidth',2);
set(gca,'FontSize',30);
%% 绘制某个刺激下反应神经元的zscore
onset=5;
zscore_response = cell(1,length(ROIIDResponsive));
zscore_inhibtion = cell(1,length(ROIIDinhibiton));
zscore_noresponse=cell(1,roi_no-length(ROIIDResponsive)-length(ROIIDinhibiton));
r=1;
nr=1;
ir=1;
for i=1:roi_no
    if ismember(i, ROIIDResponsive)==1
        zscore_response{r}=ca_zscore{i};
        r=r+1;
    elseif ismember(i,ROIIDinhibiton)==1
        zscore_inhibtion{ir}=ca_zscore{i};
        ir=ir+1;
    else
        zscore_noresponse{nr}=ca_zscore{i};
        nr=nr+1;
    end
end

ca_zscore_ppp=[];
for i=1:size(zscore_response,2) 
    dff=zscore_response{i};
    ca_zscore_ppp=cat(1,ca_zscore_ppp,dff);       
end
dff_response_mean=mean(ca_zscore_ppp);
dff_response_sem=std(ca_zscore_ppp)/sqrt(size(ca_zscore_ppp,1));

noresponse=[];
for i=1:size(zscore_noresponse,2) 
    dff=zscore_noresponse{i};
    noresponse=cat(1,noresponse,dff);       
end
noresponse_mean=mean(noresponse);
noresponse_sem=std(noresponse)/sqrt(size(noresponse,1));

inhibit_trace=[];
for i=1:size(zscore_inhibtion,2) 
    dff=zscore_inhibtion{i};
    inhibit_trace=cat(1,inhibit_trace,dff);       
end
inhibit_mean=mean(inhibit_trace);
inhibit_sem=std(inhibit_trace)/sqrt(size(inhibit_trace,1));

all_trace=[];
for i=1:roi_no
    dff=ca_zscore{i};
    all_trace=cat(1,all_trace,dff);       
end
all_trace_mean=mean(all_trace);
all_trace_sem=std(all_trace)/sqrt(size(all_trace,1));



figure;
subplot(2,4,4)
x=1:size(dff_response_ppp,2);
shadedErrorBar(x,dff_response_mean,dff_response_sem,{'-r','linewidth',2});
hold on
shadedErrorBar(x,noresponse_mean,noresponse_sem,{'-k','linewidth',2});
hold on
shadedErrorBar(x,inhibit_mean,inhibit_sem,{'-y','linewidth',2});

% figure;
% subplot(2,4,2)
% shadedErrorBar(x,all_trace_mean,all_trace_sem,{'-y','linewidth',2});

plot([onset*Fs+1,onset*Fs+1],[-2,3.5],'k--','LineWidth',2);
set(gca,'XTickLabel',{'-5','0','10'},'XTick',[1 onset*Fs+1 size(dff_response_ppp,2)]);
set(gca,'tickdir','out');     
xlim([0,size(dff_response_ppp,2)]);
ylim([-2,5.5]);
xlabel('Time from  action onset (s)','fontsize',30);
ylabel(['\Delta' 'F/F0 (%)'],'fontsize',30);
box off
set(gca,'FontSize',30);
set(gca,'LineWidth',2);
%%
zscore_ca_mean=cell(1,roi_no);
peak_value_decrease=[];
peak_value_increase=[];
mean_zscore_F_increase=[];
mean_zscore_F_decrease=[];
mean_zscore_F=cell(1,roi_no);
mean_zscore_F_b=cell(1,roi_no);
responseWindow = (floor((pre_onset)*Fs)+1):floor((pre_onset+2)*Fs);
baseWindow=(floor((base)*Fs)+1):floor((base+2)*Fs);
for i=1:roi_no
    zscore_ca_mean{i}=mean(ca_zscore{i});
%     peak_value{i}=max(zscore_ca_mean{i}(responseWindow));
    mean_zscore_F{i}=mean(zscore_ca_mean{i}(responseWindow));
    mean_zscore_F_b{i}=mean(zscore_ca_mean{i}(baseWindow));
    if ismember(i,ROIIDResponsive) 
        peak_value_increase=[peak_value_increase,max(zscore_ca_mean{i}(responseWindow))];
        mean_zscore_F_increase=[mean_zscore_F_increase,mean_zscore_F{i}];
    elseif ismember(i,ROIIDinhibiton)
        peak_value_decrease=[peak_value_decrease,min(zscore_ca_mean{i}(responseWindow))];
        mean_zscore_F_decrease=[mean_zscore_F_decrease,mean_zscore_F{i}];
    end
end

%%
%tail pinch    实验中将反应的神经元在1-neuronnum 的矩阵中标出 （手动组合后得到tail pinch_differenttrial文件）
original_matrix = ROIIDResponsive; % 原始矩阵
desired_length = size(neuronID,1); % 目标长度

neuron_matrix = zeros(1, desired_length); % 创建一个包含目标长度的0矩阵


% 将原始矩阵的值插入到指定间隔的位置

for i = 1:numel(original_matrix)
    neuron_matrix(original_matrix(i)) = original_matrix(i);
end

neuron_matrix=neuron_matrix'; % 显示补齐后的矩阵
%%  挑选指定的神经元进行zscore 绘制（目前基于tail pinch 的反应写的程序）
neuronID_response=xlsread('E:\data\Ca_analyzing\Calb1_all_ca\final\calb_responseforlist.xlsx'); %把不同的反应神经元合成了一个excel 文件再扔进来
% responseID=neuronID_response(:,1)~=0 & all(neuronID_response(:, [2:6]) == 0, 2); %clip activated
% responseID=neuronID_responseD(:,2)~=0 & all(neuronID_response(:, [1,3:6]) == 0, 2); %clip inhibtied
% responseID=neuronID_response(:,1)~=0 & neuronID_response(:,4)~=0 & all(neuronID_response(:,[2,3,5,6]) == 0, 2);  %clip activated attack inhibited
% responseID=neuronID_response(:,4)~=0 & all(neuronID_response(:,[1:3,5,6]) == 0, 2);% only attack inhibited
responseID=neuronID_response(:,5)~=0 & all(neuronID_response(:,[1:4,6]) == 0, 2);% licking activated
% responseID=neuronID_response(:,1)~=0 & neuronID_response(:,3)~=0 & neuronID_response(:,5)~=0 &all(neuronID_responseD(:,[2,4,6]) == 0, 2);
% responseID=neuronID_response(:,1)~=0 &all(neuronID_response(:,[2,3,4,5,6]) == 0, 2);
% responseID=neuronID_response(:,3)~=0 &all(neuronID_response(:,[1,2,4,5,6]) == 0, 2);
% responseID=neuronID_response(:,5)~=0 &all(neuronID_response(:,[1,2,3,4,6]) == 0, 2);
% 
responseneuron=neuronID_response(responseID,:);
non_zero_cols=any(responseneuron~=0,1);
responseneuron=responseneuron(:,non_zero_cols);
reponsneuron=responseneuron';
r=1;
for i=1:roi_no
    if ismember(i, reponsneuron)==1
        dff_designate{r}=ca_zscore{i};
        r=r+1;
    end
end
dff_designate_p=[];
for i=1:size(dff_designate,2) 
    dff=dff_designate{i};
    dff_designate_p=cat(1,dff_designate_p,dff);       
end
%% 绘制选定神经元的zscore 曲线图及热图
r=1;
for i=1:roi_no
    if ismember(i, responseID)==1
        dff_pick{r}=ca_zscore{i};
        r=r+1;
    end
end
dff_pick_p=[];
for i=1:size(dff_pick,2) 
    dff=dff_pick{i};
    dff_pick_p=cat(1,dff_pick_p,dff);       
end
dff_pick_mean=mean(dff_pick_p);
dff_pick_sem=std(dff_pick_p)/sqrt(size(dff_pick_p,1));

figure;
subplot(2,4,2)
yMax=max(dff_pick_mean + dff_pick_sem); yMin=min(dff_pick_mean - dff_pick_sem);
x=1:size(dff_transient_p,2);
shadedErrorBar(x,dff_pick_mean,dff_pick_sem,{'-r','linewidth',2});
hold on
plot([onset*Fs+1,onset*Fs+1],[-2,5.5],'k--','LineWidth',2);
set(gca,'XTickLabel',{'-5','0','10'},'XTick',[1 onset*Fs+1 size(dff_sustain_p,2)]);
set(gca,'tickdir','out');     
xlim([0,size(dff_sustain_p,2)]);
ylim([yMin,yMax]);
xlabel('Time from action onset (s)','fontsize',30);
ylabel([ 'zscore'],'fontsize',30);
box off
set(gca,'FontSize',30);
set(gca,'LineWidth',2);

%% 绘制sustain 和transient 曲线差异
r=1;
for i=1:roi_no
    if ismember(i, sustained_neun)==1
        dff_sustain{r}=ca_zscore{i};
        r=r+1;
    end
end
dff_sustain_p=[];
for i=1:size(dff_sustain,2) 
    dff=dff_sustain{i};
    dff_sustain_p=cat(1,dff_sustain_p,dff);       
end
dff_sustain_mean=mean(dff_sustain_p);
dff_sustain_sem=std(dff_sustain_p)/sqrt(size(dff_sustain_p,1));

r=1;
transient_neun=setdiff(ROIIDResponsive, sustained_neun);
for i=1:roi_no
    if ismember(i, transient_neun)==1
        dff_transient{r}=ca_zscore{i};
        r=r+1;
    end
end
dff_transient_p=[];
for i=1:size(dff_transient,2) 
    dff=dff_transient{i};
    dff_transient_p=cat(1,dff_transient_p,dff);       
end
dff_transient_mean=mean(dff_transient_p);
dff_transient_sem=std(dff_transient_p)/sqrt(size(dff_transient_p,1));

figure;
subplot(2,4,2)
yMax=max(dff_sustain_mean + dff_sustain_sem); yMin=min(dff_sustain_mean - dff_sustain_sem);
x=1:size(dff_transient_p,2);
shadedErrorBar(x,dff_sustain_mean,dff_sustain_sem,{'-r','linewidth',2});
hold on
shadedErrorBar(x,dff_transient_mean,dff_transient_sem,{'-b','linewidth',2});
hold on
plot([onset*Fs+1,onset*Fs+1],[-2,5.5],'k--','LineWidth',2);
set(gca,'XTickLabel',{'-5','0','10'},'XTick',[1 onset*Fs+1 size(dff_sustain_p,2)]);
set(gca,'tickdir','out');     
xlim([0,size(dff_sustain_p,2)]);
ylim([-2,13]);
xlabel('Time from action onset (s)','fontsize',30);
ylabel([ 'zscore'],'fontsize',30);
box off
set(gca,'FontSize',30);
set(gca,'LineWidth',2);
%% tail pinch 获得在各个情况下反应状态，0 无反应； -1 被抑制； 1 被激活
combination=unique(neuronID_response,'rows');
num_cols=size(combination,2);
for i=1:2:num_cols
     col=combination(:,i);
     no_zero_indices=col~=0;
     col(no_zero_indices)=1;
     combination(:,i)=col;
end
for i=2:2:num_cols
     col=combination(:,i);
     no_zero_indices=col~=0;
     col(no_zero_indices)=-1;
     combination(:,i)=col;
end
for i = 2:2:size(combination, 2)
    even_col = combination(:, i);
    non_zero_indices = even_col ~= 0;
    
    % 向相邻前面的奇数列同行平移
    if i > 1
        odd_col = combination(:, i - 1);
        combination(non_zero_indices, i - 1) = even_col(non_zero_indices);
        combination(non_zero_indices, i) = 0;
    end
end
combination(:,2:2:6)=[];

% 获取矩阵中不同行的种类和每种行出现的次数
[unique_rows, ~, idx] = unique(combination, 'rows');
counts = accumarray(idx, 1);

% 按照行出现的次数降序排列
[sorted_counts, sorted_indices] = sort(counts, 'descend');
sorted_unique_rows = unique_rows(sorted_indices, :);

% 构建新的矩阵
combination_p = [];
for i = 1:length(sorted_counts)
    combination_p = [combination_p; repmat(sorted_unique_rows(i, :), sorted_counts(i), 1)];
end
