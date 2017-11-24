function  LRAs=Improved_LRA( pointcloud, indices, Support_radius)
%  input: pointcloud - key_point cloud.
%         indices - the indices of key points on the key_point cloud.
%         Support_radius - the support radius for determining the local surface to calculate LRAs. 
%
%  output: LRAs - the LRAs generated on the key points.
%
%  Author: Bao Zhao {zhaobao1988@sjtu.edu.cn}

L=length(indices);
LRAs=zeros(L,3);
Index=rangesearch(pointcloud,pointcloud(indices,:),Support_radius);
for i=1:L
    
    key_point=pointcloud(indices(i),:);   
    local_points=pointcloud(Index{i,1},:); 
    if length(Index{i,1})<10  % If the radius neighbors is less than 10 points, the LRAs is directly assigned with [0,0,1]. 
        LRAs(i,:)=[0,0,1];
    end
    
    %============================obtain the direction of LRAs============================%
    mean_point=mean(local_points);
    ml_pointcloud=[local_points(:,1)-mean_point(1), local_points(:,2)-mean_point(2), local_points(:,3)-mean_point(3)];
    Cov=ml_pointcloud'*ml_pointcloud/length(local_points);
    [V,D]=eig(Cov);
    eigenvalue=[D(1,1),D(2,2),D(3,3)];
    min_d=(eigenvalue==min(eigenvalue));
    Z=V(:,min_d);
    
    %============================disambiguating the sign of LRAs============================%
    mmm=[local_points(:,1)-key_point(1), local_points(:,2)-key_point(2), local_points(:,3)-key_point(3)];
    sum_local=sum(Z'*mmm');
       if sum_local<0 
           Z=-Z;
       end
    LRAs(i,:)=Z;
    
end
end

