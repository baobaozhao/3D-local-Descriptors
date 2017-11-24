function index=judge_circle_location(circle_step,dis )
k=length(circle_step);
% if dis<circle_step(1)
% 
%     index=1;
%     return;
% 
% end
for i=1:k
    if dis<=circle_step(i)
        index=i;
        return
    end
end
index=k;
end

