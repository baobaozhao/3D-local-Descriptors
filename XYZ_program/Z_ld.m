 function Hist_out=Z_ld(pointcloud,Key_indices,LPAs,R,Z_size,deviation_angle_size,LRFs)
%  Author: Bao Zhao {zhaobao1988@sjtu.edu.cn}
%  input: pointcloud - key_point cloud.
%         Key_indices - the neighbor_indices of key points on the key_point cloud.
%         LPAs - LPA of the pointcloud.
%         R - the support radius for determining the local surface.
%         Z_size - the dimension of Z direction.
%         deviation_angle_size - the dimension of deviation angle between LPA and LRA.
%         LRAs - the LRA Attached to each key key_point.

%  output:Hist_out - the vectors of Z_ld descriptor generated on all key pointcloud.
L_key_points=length(Key_indices);
hist_size=deviation_angle_size*Z_size;
Hist=zeros(L_key_points,hist_size);
deviation_angle_step=180/deviation_angle_size;
Z_step=2*R/Z_size;
ad_index=rangesearch(pointcloud,pointcloud(Key_indices,:),R);


for i=1:L_key_points
    hist=zeros(1,hist_size);
    LRA=LRFs(i,7:9);
    key_point=pointcloud(Key_indices(i),:);
    neighbor_indices=ad_index{i,1};
    
    neighbor_indices=neighbor_indices(2:end);  % removing the index of the key point.
    neigh_size=length(neighbor_indices);
    
    if (neigh_size<10)  % if radius neighbors are less than 10 points, we randomly generate a descriptor vector at that key point.
        Hist(i,:)=rand(1,hist_size);
        Hist(i,:)=Hist(i,:)/sum(Hist(i,:));
        continue;
    end
    
    local_points=pointcloud(neighbor_indices,:);
    local_points=[local_points(:,1)-key_point(1),local_points(:,2)-key_point(2),...
        local_points(:,3)-key_point(3)];
    
    % transforming local points.
    M1=[LRFs(i,1:3);LRFs(i,4:6);LRFs(i,7:9)];
    local_points=(M1*local_points')';

    
    Cos_zhi=zeros(neigh_size,1);
    
    %Normal_length=zeros(neigh_size,1);
    
    for j=1:neigh_size
        Cos_zhi(j)=LPAs(neighbor_indices(j),:)*LRA';
        
    %Normal_length(j)=local_points(j,:)*LRA';
        
        if Cos_zhi(j)>=0.99999
            Cos_zhi(j)=0;
        elseif Cos_zhi(j)<=-0.9999
            Cos_zhi(j)=180;
        else
           Cos_zhi(j)=180*acos(Cos_zhi(j))/pi;
        end
        
        deviation_angle_list=ceil((Cos_zhi(j))/deviation_angle_step);
        
        Z_list=ceil((local_points(j,3)+R)/Z_step);
        
        if deviation_angle_list==0
            deviation_angle_list=1;
        end
        if deviation_angle_list>deviation_angle_size
            deviation_angle_list=deviation_angle_size;
        end
        
        if Z_list==0
            Z_list=1;
        end
        if Z_list>Z_size
            Z_list=Z_size;
        end
        
        index=(Z_list-1)*deviation_angle_size+deviation_angle_list;
        hist(index)=hist(index)+1;
        
    end
    Hist(i,:)=hist/neigh_size;
    
end
Hist_out=Hist;
 end