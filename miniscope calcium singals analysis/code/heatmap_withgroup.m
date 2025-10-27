% load('E:\data\tnc miniscope\Tnc-all\final\dffstructtnctail.mat');
% get_dff=dff_data;
% attack_activation=load('E:\data\tnc miniscope\Tnc-all\final\attackactivated.mat');
% attack_activation=attack_activation.ROIIDResponsive;
% attack_inhibition=load('E:\data\tnc miniscope\Tnc-all\final\attackinhibited.mat');
% attack_inhibition=attack_inhibition.ROIIDinhibiton;
% clip_activationID=load('E:\data\tnc miniscope\Tnc-all\final\pinchactivated.mat');
% clip_activationID=clip_activationID.ROIIDResponsive;
% clip_inhibitionID=load('E:\data\tnc miniscope\Tnc-all\final\pinchinhibited.mat');
% clip_inhibitionID=clip_inhibitionID.ROIIDinhibiton;
% licking_activationID=load('E:\data\tnc miniscope\Tnc-all\final\lickactivated.mat');
% licking_activationID=licking_activationID.ROIIDResponsive;
% licking_inhibition=load('E:\data\tnc miniscope\Tnc-all\final\lickinhibited.mat');
% licking_inhibition=licking_inhibition.ROIIDinhibiton;
% [tail_licking,tail_l_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'l','2');
% [tail_attacking,tail_a_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'u','2');%,'2');
% [tail_clip,tail_j_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'j','2');
%%
load('E:\data\Ca_analyzing\Calb1_all_ca\final\dffstructcalb1tail.mat');
get_dff=dff_data;
attack_activation=load('E:\data\Ca_analyzing\Calb1_all_ca\final\attackactivated.mat');
attack_activation=attack_activation.ROIIDResponsive;
attack_inhibition=load('E:\data\Ca_analyzing\Calb1_all_ca\final\attackinhibited.mat');
attack_inhibition=attack_inhibition.ROIIDinhibiton;
clip_activationID=load('E:\data\Ca_analyzing\Calb1_all_ca\final\pinchactivated.mat');
clip_activationID=clip_activationID.ROIIDResponsive;
clip_inhibitionID=load('E:\data\Ca_analyzing\Calb1_all_ca\final\pinchinhibited.mat');
clip_inhibitionID=clip_inhibitionID.ROIIDinhibiton;
licking_activationID=load('E:\data\Ca_analyzing\Calb1_all_ca\final\lickactivated.mat');
licking_activationID=licking_activationID.ROIIDResponsive;
licking_inhibition=load('E:\data\Ca_analyzing\Calb1_all_ca\final\lickinhibited.mat');
licking_inhibition=licking_inhibition.ROIIDinhibiton;
[tail_licking,tail_l_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'l','2');
[tail_attacking,tail_a_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'u','2');%,'2');
[tail_clip,tail_j_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'j','2');
%%
% tail_attacking_all=[tail_attacking_tnc,tail_attacking_calb];
% tail_licking_all=[tail_licking_tnc,tail_licking_calb];
%%
%tnc 读取后执行
% roi_no=size(tail_licking,2);
% lick_response_tnc = cell(1,length(licking_activationID));
% attack_response_tnc=cell(1,length(attack_activation));
% lr=1;
% ar=1;
% for i=1:roi_no
%     if ismember(i, licking_activationID)==1
%         lick_response_tnc{lr}=tail_licking{i};
%         lr=lr+1;
%     elseif ismember(i,attack_activation)==1
%         attack_response_tnc{ar}=tail_attacking{i};
%         ar=ar+1;
%     end
% end

%calb 读取后执行
% roi_no=size(tail_licking,2);
% lick_response_calb = cell(1,length(licking_activationID));
% attack_response_calb=cell(1,length(attack_activation));
% lr=1;
% ar=1;
% for i=1:roi_no
%     if ismember(i, licking_activationID)==1
%         lick_response_calb{lr}=tail_licking{i};
%         lr=lr+1;
%     elseif ismember(i,attack_activation)==1
%         attack_response_calb{ar}=tail_attacking{i};
%         ar=ar+1;
%     end
% end
%%
% 上面执行后合并
% lick_response=[lick_response_tnc,lick_response_calb];%Tnc 在前 calb 在后
% attack_response=[attack_response_tnc,attack_response_calb];%Tnc 在前 calb 在后
% lick_average=[];
% for i=1:size(lick_response,2) 
%     dff=lick_response{i};
%     lick_average=cat(1,lick_average,dff);       
% end
% lick_average_mean=mean(lick_average);
% attack_average=[];
% for i=1:size(attack_response,2) 
%     dff=attack_response{i};
%     attack_average=cat(1,attack_average,dff);       
% end
% attack_average_mean=mean(attack_average);
%%
pre_onset=5;%反应起点
Fs = 14.854;
baseline=int32(5*Fs);
zscore_lick=calculate_zscore(tail_licking,baseline);
zscore_attack=calculate_zscore(tail_attacking,baseline);
zscore_clip=calculate_zscore(tail_clip,baseline);
noresponse=[];
allactive=[];
allinhibitive=[];

numberofneun=[];

for i=1:length(tail_l_rawID)   
    if ismember(i,clip_activationID) && ismember(i,licking_activationID)&& ismember(i,attack_activation)
        allactive=[allactive,i];
    elseif ismember(i,clip_inhibitionID) && ismember(i,licking_inhibition)&& ismember(i,attack_inhibition)
        allinhibitive=[allinhibitive,i];
    end
end
inum=[allactive,allinhibitive];

roi=1:length(tail_l_rawID);
roi=setdiff(roi,inum);
indexcl_a=pick_neuron(roi,clip_activationID,licking_activationID,1);
roi=setdiff(roi,indexcl_a);
indexcl_i=pick_neuron(roi,clip_inhibitionID,licking_inhibition,1);
roi=setdiff(roi,indexcl_i);
inum=[inum,indexcl_a,indexcl_i];
indexca_a=pick_neuron(roi,clip_activationID,attack_activation,1);
roi=setdiff(roi,indexca_a);
indexca_i=pick_neuron(roi,clip_inhibitionID,attack_inhibition,1);
roi=setdiff(roi,indexca_i);
inum=[inum,indexca_a,indexca_i];
indexcl_ocla=pick_neuron(roi,clip_activationID,licking_inhibition,1);
roi=setdiff(roi,indexcl_ocla);
indexcl_ocli=pick_neuron(roi,clip_inhibitionID,licking_activationID,1);
roi=setdiff(roi,indexcl_ocli);
inum=[inum,indexcl_ocla,indexcl_ocli];
indexca_ocaa=pick_neuron(roi,clip_activationID,attack_inhibition,1);
roi=setdiff(roi,indexca_ocaa);
indexca_ocai=pick_neuron(roi,clip_inhibitionID,attack_activation,1);
roi=setdiff(roi,indexca_ocai);
inum=[inum,indexca_ocaa,indexca_ocai];


pinch_only=pick_neuron(roi,clip_activationID,inum,0);
pinch_only=[pinch_only,pick_neuron(roi,clip_inhibitionID,inum,0)];
inum=[pinch_only,inum];
attack_only=pick_neuron(roi,attack_activation,inum,0);
attack_only=[attack_only,pick_neuron(roi,attack_inhibition,inum,0)];
inum=[attack_only,inum];

lick_only=pick_neuron(roi,licking_activationID,inum,0);
lick_only=[lick_only,pick_neuron(roi,licking_inhibition,inum,0)];
inum=[lick_only,inum];

noresponse=setdiff(roi,inum);
index=[lick_only,indexcl_a,indexcl_i,indexcl_ocli,indexcl_ocla,attack_only,indexca_a,indexca_i,indexca_ocai,indexca_ocaa,pinch_only,allactive,allinhibitive,noresponse];
index=flip(index);
numberofneun=[numberofneun,length(lick_only),length(indexcl_a)+length(indexcl_i),length(indexcl_ocla)+length(indexcl_ocli)];
numberofneun=[numberofneun,length(attack_only),length(indexca_a)+length(indexca_i),length(indexca_ocaa)+length(indexca_ocai)];
numberofneun=[numberofneun,length(pinch_only),length(allactive)+length(allinhibitive),length(noresponse)];
%%
dff_analyzing=tail_clip;
ca_zscore=calculate_zscore(dff_analyzing,baseline);
onset=5;
%%
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

zscore_plot_sorted = zscore_plot(index,:);
% 使用imagesc函数绘制热图
figure('Position', [100, 100, 400, 800]); 
imagesc(zscore_plot_sorted, [N,M]);
hold on
% 
for i = 1:length(numberofneun)
    % 在每个 y 坐标处绘制一条横线
    y =length(tail_l_rawID)-sum(numberofneun(1:i))+0.5;
    plot([0, size(zscore_plot,2)], [y, y], 'k-', 'LineWidth', 1.5); % 横线
end

plot([onset*Fs+1,onset*Fs+1],[0,size(zscore_plot,1)+0.5],'y--','LineWidth',2);
set(gca,'XTickLabel',{'-5','0','10'},...
    'XTick',[1 onset*Fs+1 size(zscore_plot,2)]);
set(gca,'tickdir','out');  
% set(gca,'Ytick',[],'YtickLabel',[]);
xlabel('Time from event onset (s)','fontsize',30);
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
caxis([-2,2]);
actucolorbar = colorbar;
% set(actucolorbar,'position',[0.51 0.7 0.006 0.15]);
title(actucolorbar,['  zscore'],'fontsize',25);
set(gca,'LineWidth',2);
set(gca,'FontSize',30);
%% zscore 函数
function ca_zscore=calculate_zscore(dff_analyzing,baseline)
        cell_means = cellfun(@(x) mean(x(:, 1:baseline), 2), dff_analyzing, 'UniformOutput', false);
        cell_stds = cellfun(@(x) std(x(:, 1:baseline), 0, 2), dff_analyzing, 'UniformOutput', false);
        ca_zscore = cellfun(@(x, m, s) (x - m) ./ s, dff_analyzing, cell_means, cell_stds, 'UniformOutput', false);       
end
function neuronID=pick_neuron(testlist,neuronlistA,neuronlistB,choice)
     neuronID=[];
     if choice==1
         for i=1:length(testlist)
             roi=testlist(i);
             if ismember(roi,neuronlistA) && ismember(roi,neuronlistB)
                 neuronID=[neuronID,roi];
             end
         end
     elseif choice==0
         for i=1:length(testlist)
             roi=testlist(i);
             if ismember(roi,neuronlistA) && ~ismember(roi,neuronlistB)
                 neuronID=[neuronID,roi];
             end
         end
     end
end
