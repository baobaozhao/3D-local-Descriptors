function Hist_out=XYZ_spatial(pointcloud,Key_indices,R,X_size,Y_size,Z_size,LRFs)
%  Author: Bao Zhao {zhaobao1988@sjtu.edu.cn}
%  input: pointcloud - key_point cloud.
%         Key_indices - the neighbor_indices of key points on the key_point cloud.
%         R - the support radius for determining the local surface.
%         X_size - the dimension along X direction.
%         Y_size - the dimension along Y direction.
%         Z_size - the dimension along Z direction.
%         LRFs - the LRF Attached to each key key_point.

%  output:Hist_out - the vectors of XYZ_spatial descriptor generated on all key pointcloud.

L_key_points=length(Key_indices);
hist_size=X_size*Y_size*Z_size;

Hist=zeros(L_key_points,hist_size);

X_step=2*R/X_size;
Y_step=2*R/Y_size;
Z_step=2*R/Z_size;

ad_index=rangesearch(pointcloud,pointcloud(Key_indices,:),R);

for i=1:L_key_points
    hist=zeros(1,hist_size);
    %LRA=LRFs(i,7:9);
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
    
    
    for j=1:neigh_size
        
        
        X_list=ceil((local_points(j,1)+R)/X_step);
        Y_list=ceil((local_points(j,2)+R)/Y_step);
        Z_list=ceil((local_points(j,3)+R)/Z_step);
        
        
        if X_list==0
            X_list=1;
        end
        if X_list>X_size
            X_list=X_size;
        end
        
        if Y_list==0
            Y_list=1;
        end
        if Y_list>Y_size
            Y_list=Y_size;
        end
        
        if Z_list==0
            Z_list=1;
        end
        if Z_list>Z_size
            Z_list=Z_size;
        end
        
        index=(X_list-1)*Y_size*Z_size+(Y_list-1)*Z_size+Z_list;
        
        hist(index)=hist(index)+1;
        
    end
    Hist(i,:)=hist/neigh_size;
    
end
Hist_out=Hist;
 end