load('E:\data\Ca_analyzing\Calb1_all_ca\final\dffstructcalb1tail.mat')
sample_no=length (unique (fieldnames (dffstruct)));
get_dff=dff_data;
[roi_no]=get_dff.getroi_no(dffstruct);
responseneuron=zeros(sample_no,5);
for i=1:sample_no
    responseneuron(i,1)=i;  
end
%根据行为筛选
% [tail_licking,tail_l_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'l','2');
% [tail_attacking,tail_a_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'u','2')%,'2');
[tail_clip,tail_j_rawID]=get_dff.getneuronlistbytrialtype(dffstruct,'j','2');

%%
dff_analyzing=tail_clip; %挑选分析的数据
rawID=tail_j_rawID;
%%
Fs = 14.854; 
onset=5;
baseline=int32(5*Fs);
cell_means = cellfun(@(x) mean(x(:, 1:baseline), 2), dff_analyzing, 'UniformOutput', false);
cell_stds = cellfun(@(x) std(x(:, 1:baseline), 0, 2), dff_analyzing, 'UniformOutput', false);

ca_zscore = cellfun(@(x, m, s) (x - m) ./ s, dff_analyzing, cell_means, cell_stds, 'UniformOutput', false);
%%  挑选指定的神经元进行zscore 绘制（目前基于tail pinch 的反应写的程序）
neuronID=xlsread('E:\data\Ca_analyzing\Calb1_all_ca\final\calb_responseforlist.xlsx'); %把不同的反应神经元合成了一个excel 文件再扔进来
% responseID=neuronID(:,1)~=0 & all(neuronID(:, [3,5]) == 0, 2); %clip activated
% responseID=neuronID(:,2)~=0 & all(neuronID(:, [1,3:6]) == 0, 2); %clip inhibtied
% responseID=neuronID(:,1)~=0 & neuronID(:,4)~=0 & all(neuronID(:,[2,3,5,6]) == 0, 2);  %clip activated attack inhibited
% responseID=neuronID(:,3)~=0 & all(neuronID(:,[1,5]) == 0, 2);% attack activated
% responseID=neuronID(:,5)~=0 & all(neuronID(:,[1,3]) == 0, 2);% licking activated
% responseID=neuronID(:,1)~=0 & neuronID(:,3)~=0 & neuronID(:,5)~=0 &all(neuronID(:,[2,4,6]) == 0, 2);
% responseID=neuronID(:,1)~=0 &all(neuronID(:,[2,3,4,5,6]) == 0, 2);
% responseID=neuronID(:,3)~=0 &all(neuronID(:,[1,2,4,5,6]) == 0, 2);
% responseID=neuronID(:,5)~=0 &all(neuronID(:,[1,2,3,4,6]) == 0, 2);
% responseID=neuronID(:,1)~=0 & neuronID(:,5)~=0 & all(neuronID(:,[3]) == 0, 2);  %clip activated lick activated
% responseID=neuronID(:,1)~=0 & neuronID(:,3)~=0 & all(neuronID(:,[5]) == 0, 2);  %clip activated attack activated
responseID=neuronID(:,1)~=0 & neuronID(:,3)~=0 & neuronID(:,5)~=0;
responseneuron=neuronID(responseID,:);
non_zero_cols=any(responseneuron~=0,1);
responseneuron=responseneuron(:,non_zero_cols);
reponsneuron=responseneuron';

%% 绘制选定神经元的zscore 曲线图及热图
r=1;
for i=1:roi_no
    if ismember(i, reponsneuron)==1
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
% yMax=max(dff_pick_mean + dff_pick_sem); yMin=min(dff_pick_mean - dff_pick_sem);
x=1:size(dff_pick_p,2);
shadedErrorBar(x,dff_pick_mean,dff_pick_sem,{'-r','linewidth',2});
hold on
plot([onset*Fs+1,onset*Fs+1],[-2,5.5],'k--','LineWidth',2);
set(gca,'XTickLabel',{'-5','0','10'},'XTick',[1 onset*Fs+1 size(dff_pick_p,2)]);
set(gca,'tickdir','out');     
xlim([0,size(dff_pick_p,2)]);
ylim([-0.5,6]);
xlabel('Time from action onset (s)','fontsize',30);
ylabel([ 'zscore'],'fontsize',30);
box off
set(gca,'FontSize',30);
set(gca,'LineWidth',2);