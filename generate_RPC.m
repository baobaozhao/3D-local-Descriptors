function [Pre,Recall]=generate_RPC(Model,Scene,Model_Hist,Scene_Hist,support_radius)

[N,d]=knnsearch(Scene_Hist,Model_Hist,'k',2);
Ratio_hist=d(:,1)./d(:,2);

[Num,~]=size(Scene_Hist);
[~,In]=sort(Ratio_hist(:,1));
Pre=zeros(10,1);
Recall=zeros(10,1);
k=1;
for i=0.1:1.1:10
     num=round(length(Ratio_hist)*(i*0.1));
     Whole_index=N(In(1:num),1);
     Part_index=In(1:num);
     Dis=((Model(Part_index,1)-Scene(Whole_index,1)).^2+(Model(Part_index,2)-Scene(Whole_index,2)).^2 ...
         +(Model(Part_index,3)-Scene(Whole_index,3)).^2).^(0.5);
     list=find(Dis<(support_radius/2));
    Pre(k)=(length(Whole_index)-length(list))/length(Whole_index);
    Recall(k)=(length(list))/Num;
    k=k+1;
end

end