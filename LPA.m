function LPAs= LPA(pointcloud, indices, r)

%  input: pointcloud - point cloud.
%         indices - the indices of key pointcloud on the point cloud.
%         r - the support radius for determining the local surface to calculate LRA. 
% 
%  output:LPAs - the LPA generated on the key points.
%
%  Author: Bao Zhao {zhaobao1988@sjtu.edu.cn}

L=length(indices);
ad_index=rangesearch(pointcloud,pointcloud(indices,:),r);
LPAs=zeros(L,3);

for i=1:L
    point=pointcloud(indices(i),:);   
    
    if (length(ad_index{i,1})<10) % if the radius neighbors are less than 10 points, the 10 closest points are used.  
         [N,~]=knnsearch(pointcloud,pointcloud(i,:),'k',10);
         local_pointcloud=pointcloud(N,:);
    else
        local_pointcloud=pointcloud(ad_index{i,1},:); 
    end
    
    
    
    %============================obtaining the direction of LPA============================%
    mean_point=mean(local_pointcloud);
    ml_pointcloud=[local_pointcloud(:,1)-mean_point(1), local_pointcloud(:,2)-mean_point(2), local_pointcloud(:,3)-mean_point(3)];
    Cov=ml_pointcloud'*ml_pointcloud/length(local_pointcloud);
    [V,D]=eig(Cov);
    eigenvalue=[D(1,1),D(2,2),D(3,3)];
    min_d=(eigenvalue==min(eigenvalue));
    Z=V(:,min_d);
    
    
    %============================disambiguating the sign of LPA============================%
    mmm=[local_pointcloud(:,1)-point(1), local_pointcloud(:,2)-point(2), local_pointcloud(:,3)-point(3)];
    sum_local=sum(Z'*mmm');
    if sum_local<0 
        Z=-Z;
    end
    LPAs(i,:)=Z/norm(Z);
    
end
end

