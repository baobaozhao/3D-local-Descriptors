 function Hist_out=radial_ld(pointcloud,Key_indices,LPAs,R,radial_size,deviation_angle_size,LRAs)
%  Author: Bao Zhao {zhaobao1988@sjtu.edu.cn}
%  input: pointcloud - key_point cloud.
%         Key_indices - the neighbor_indices of key points on the key_point cloud.
%         LPAs - LPA of the pointcloud.
%         R - the support radius for determining the local surface.
%         radial_size - the dimension along radial direction.
%         deviation_angle_size - the dimension of deviation angle between LPA and LRA.
%         LRAs - the LRA Attached to each key key_point.

%  output:Hist_out - the vectors of radial_ld descriptor generated on all key pointcloud.

L_key_points=length(Key_indices);
hist_size=deviation_angle_size*radial_size;


Hist=zeros(L_key_points,hist_size);
deviation_angle_step=180/deviation_angle_size;
radial_step=R/radial_size;

ad_index=rangesearch(pointcloud,pointcloud(Key_indices,:),R);

for i=1:L_key_points
    hist=zeros(1,hist_size);
    LRA=LRAs(i,:);
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
%     dirct=cross(LRA,[0,0,1]); % the direction of transforming.
%     theta=acos(LRA*[0;0;1]); % the angle of transforming.
%     dirct=dirct(:)/norm(dirct);
%     cosine=cos(theta);
%     sine=sin(theta);
%     x=dirct(1);
%     y=dirct(2);
%     z=dirct(3);
%     M1=[cosine+(1-cosine)*x*x,   (1-cosine)*x*y-sine*z,    (1-cosine)*x*z+sine*y
%         (1-cosine)*y*x+sine*z,   cosine+(1-cosine)*y*y,    (1-cosine)*y*z-sine*x
%         (1-cosine)*z*x-sine*y,   (1-cosine)*z*y+sine*x,  cosine+(1-cosine)*z*z];
%     local_points=(M1*local_points')';

    
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
        Sphere_L=sqrt(local_points(j,1)^2+local_points(j,2)^2+local_points(j,3)^2);
        
        
        %radial_step=R/radial_size;
        radial_list=ceil(Sphere_L/radial_step);
        
        
        if deviation_angle_list==0
            deviation_angle_list=1;
        end
        if deviation_angle_list>deviation_angle_size
            deviation_angle_list=deviation_angle_size;
        end
        
        if radial_list==0
            radial_list=1;
        end
        if radial_list>radial_size
            radial_list=radial_size;
        end
        
        index=(radial_list-1)*deviation_angle_size+deviation_angle_list;
        
        hist(index)=hist(index)+1;
    end
    Hist(i,:)=hist/neigh_size;
    
end

Hist_out=Hist;

 end