%%
% 通过不同条件筛选需要的神经元的dff; trialtype例:'l';'f';stimtyep:'1','2'...;genetype:'1','2'...
%%
function get_dff=dff_data
       get_dff.getneuronlistbytrialtype=@getneuronlistbytrialtype;
       get_dff.getneuronlistbygenetype=@getneuronlistbygenetype;
       get_dff.getroi_no=@getroi_no;
end
function [roi_no]=getroi_no(dffstruct)
     neuronID=fieldnames(dffstruct);    
     for i=1:size(neuronID,1) 
           mouseID{i}=str2num(neuronID{i}(4:5));   
     end
     mouseID=cell2mat(mouseID);
     [mouseID,~]=unique(mouseID,'stable');
     roi_no=0;
     for n=1:size(mouseID,2)
         roiID=cell(1,200);
         i=1;
         for m=1:size(neuronID,1)
             if neuronID{m}(4:5)==num2str(mouseID(n))
                 roiID{i}=str2num(neuronID{m}(8:length(neuronID{m})));
                 i=i+1;
             end
         end
         roiID=cell2mat(roiID);
         roi_no=roi_no+max(roiID);
     end
end
function [dfftrials,rawID]=getneuronlistbytrialtype(dffstruct,trialtype,stimtype,genetype)
    if nargin<4
        neuronID=fieldnames(dffstruct);
        rawID=[];
        j=1;
        for i= 1: size(neuronID,1)
            if neuronID{i}(1)==trialtype && neuronID{i}(6)==stimtype
               dfftrials{j}=dffstruct.(neuronID{i});
               rawID=cat(1,rawID,neuronID{i});
               j=j+1;
            end
        end
    elseif nargin<5
        neuronID=fieldnames(dffstruct);
        rawID=[];
        j=1;
        for i= 1: size(neuronID,1)
            if neuronID{i}(1)==trialtype && neuronID{i}(6)==stimtype && neuronID{i}(2)==genetype
               dfftrials{j}=dffstruct.(neuronID{i});
               rawID=cat(1,rawID,neuronID{i});
               j=j+1;
            end
        end
    end
end
function [dfftrials]=getneuronlistbygenetype(dffstruct,genetype,roi_no)
    neuronID=fieldnames(dffstruct);
    dfftrials=cell(1,roi_no);
    for j=1:roi_no
        for i= 1: size(neuronID,1)
            if neuronID{i}(2)==genetype
               dfftrials{j}=dffstruct.(neuronID{i});
            end
        end
    end
end
