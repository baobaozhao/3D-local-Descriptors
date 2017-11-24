function Hist_out=radial_azimuth_elevation(pointcloud,Key_indices,R,radial_size,azimuth_angle_size,elevation_angle_size,LRFs)
%  Author: Bao Zhao {zhaobao1988@sjtu.edu.cn}
%  input: pointcloud - key_point cloud.
%         Key_indices - the neighbor_indices of key points on the key_point cloud.
%         LPAs - LPA of the pointcloud.
%         R - the support radius for determining the local surface.
%         radial_size - the dimension along radial direction.
%         azimuth_angle_size - the dimension along azimuth direction.
%         elevation_angle_size - the dimension along elevation direction.
%         deviation_angle_size - the dimension of deviation angle between LPA and LRA.
%         LRFs - the LR Attached to each key key_point.

%  output:Hist_out - the vectors of radial_azimuth_elevation descriptor generated on all key pointcloud.

L_key_points=length(Key_indices);
hist_size=radial_size*azimuth_angle_size*elevation_angle_size;

Hist=zeros(L_key_points,hist_size);

%X_step=R/radial_size; %  球径向距离。

radial_step_even=zeros(radial_size,1);
for i=1:radial_size
    radial_step_even(i)=(i)^(1/3);
end
radial_rxi=R/radial_step_even(end);
radial_step_even=radial_step_even*radial_rxi;

Y_angle_step=180/azimuth_angle_size; %  在xy平面中的角度角度计算。

elevation_angle_step=zeros(elevation_angle_size,1);
for i=1:elevation_angle_size
    
    end_degree=1-i*(2/elevation_angle_size);
    elevation_angle_step(i)=acos(end_degree)*180/pi;
    
end

[ad_index,Dist]=rangesearch(pointcloud,pointcloud(Key_indices,:),R);

for i=1:L_key_points
    hist=zeros(1,hist_size);
    %LRA=LRFs(i,7:9);
    key_point=pointcloud(Key_indices(i),:);
    neighbor_indices=ad_index{i,1};
    
    dists=Dist{i,1};
    dists=dists(2:end);
    
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
        
        radial_list=judge_circle_location(radial_step_even,dists(j));
        
        azimuth_angle_list=find_y_angle_list(local_points(j,:),Y_angle_step);
        
        if dists(j)~=0
        z_angle=180*acos(local_points(j,3)/dists(j))/pi;
        elevation_angle_list=judge_circle_location(elevation_angle_step,z_angle);
        else
            elevation_angle_list=1;
        end
       % radial_list=ceil((local_points(j,1)+R)/X_step);
        %azimuth_angle_list=ceil((local_points(j,2)+R)/Y_step);
        %Z_list=ceil((local_points(j,3)+R)/Z_step);
        
        if radial_list==0 
            radial_list=1;
        end
        
        if radial_list>radial_size
            radial_list=radial_size;
        end
        
        if azimuth_angle_list==0
            azimuth_angle_list=1;
        end
        
        if azimuth_angle_list>azimuth_angle_size
            azimuth_angle_list=azimuth_angle_size;
        end
        
        if elevation_angle_list==0
            elevation_angle_list=1;
        end
        
        if elevation_angle_list>elevation_angle_size
            elevation_angle_list=elevation_angle_size;
        end
        
%         local_points(j,:)
%         azimuth_angle_list
%         elevation_angle_list
        
        
        index=(radial_list-1)*azimuth_angle_size*elevation_angle_size+(azimuth_angle_list-1)*elevation_angle_size+elevation_angle_list;
        
        if isnan(index)
            
            bao=1;
            
        end
        
        
        if index>hist_size||index<1
            
            bao=1;
            
        end
        
        
        
        hist(index)=hist(index)+1;
        
    end
    
    Hist(i,:)=hist/neigh_size;
    
end
Hist_out=Hist;
 end