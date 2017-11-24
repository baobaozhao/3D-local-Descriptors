 function Hist_out=SDASS(pointcloud,Key_indices,LPAs,R,height_size,projected_radial_size,deviation_angle_size,LRAs)
%  Author: Bao Zhao {zhaobao1988@sjtu.edu.cn}
%  input: pointcloud - key_point cloud.
%         Key_indices - the neighbor_indices of key points on the key_point cloud.
%         LPAs - LPA of the pointcloud.
%         R - the support radius for determining the local surface to calculate SDASS descriptor.
%         projected_radial_size - the dimension of height distance.
%         height_size - the dimension of height distance.
%         deviation_angle_size - the dimension of deviation angle between LPA and LRA.
%         LRAs - the LRA Attached to each key key_point.

%  output:Hist_out - the vectors of SDASS descriptor generated on all key pointcloud.

L_key_points=length(Key_indices);
hist_size=deviation_angle_size*height_size*projected_radial_size;

Hist=zeros(L_key_points,hist_size);
deviation_angle_step=180/deviation_angle_size;
height_step=2*R/height_size;
projected_radial_step=R/projected_radial_size;
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
     %theta=acos(LRA*[0;0;1]); % the angle of transforming.
    theta=LRA*[0;0;1];
     
     
   if theta<-0.999
        
        local_points=-local_points;
        
    elseif theta<0.999&&theta>-0.999
        theta=acos(theta);
         dirct=cross(LRA,[0,0,1]); % the direction of transforming.
        dirct=dirct(:)/norm(dirct);
        cosine=cos(theta);
        sine=sin(theta);
        x=dirct(1);
        y=dirct(2);
        z=dirct(3);
        M1=[cosine+(1-cosine)*x*x,   (1-cosine)*x*y-sine*z,    (1-cosine)*x*z+sine*y
            (1-cosine)*y*x+sine*z,   cosine+(1-cosine)*y*y,    (1-cosine)*y*z-sine*x
            (1-cosine)*z*x-sine*y,   (1-cosine)*z*y+sine*x,  cosine+(1-cosine)*z*z];
        
        local_points=(M1*local_points')';
        
    end
    

    Cos_zhi=zeros(neigh_size,1);
    Normal_length=zeros(neigh_size,1);
    for j=1:neigh_size
        Cos_zhi(j)=LPAs(neighbor_indices(j),:)*LRA';
        Normal_length(j)=local_points(j,:)*LRA';
        
        if Cos_zhi(j)>=0.99999
            Cos_zhi(j)=0;
        elseif Cos_zhi(j)<=-0.9999
            Cos_zhi(j)=180;
        else
           Cos_zhi(j)=180*acos(Cos_zhi(j))/pi;
        end
        
        deviation_angle_list=ceil((Cos_zhi(j))/deviation_angle_step);
        height_list=ceil((local_points(j,3)+R)/height_step);
        projected_radial_L=sqrt(local_points(j,1)^2+local_points(j,2)^2);
        projected_radial_list=ceil(projected_radial_L/projected_radial_step);
        
        
        if deviation_angle_list==0
            deviation_angle_list=1;
        end
        if deviation_angle_list>deviation_angle_size
            deviation_angle_list=deviation_angle_size;
        end
        if height_list==0
            height_list=1;
        end
        if height_list>height_size
            height_list=height_size;
        end
        if projected_radial_list==0
            projected_radial_list=1;
        end
        if projected_radial_list>projected_radial_size
            projected_radial_list=projected_radial_size;
        end
        
        index=(height_list-1)*projected_radial_size*deviation_angle_size+(projected_radial_list-1)*deviation_angle_size+deviation_angle_list;
        
        hist(index)=hist(index)+1;
        
    end
    Hist(i,:)=hist/neigh_size;
    Hist_out=Hist;
end

%============================eliminating redundant bins from the histogram============================%
% Num=floor(height_size/2);
% Index_eliminating=[];
% for i=1:Num
%     L_r=R-sqrt(R^2-(R-i*height_step)^2);
%     ratio_zhi=L_r/projected_radial_step;
%     if ratio_zhi<1
%         break
%     else
%         num_radial=floor(ratio_zhi);
%         for j=1:num_radial
%             Index_eliminating=[Index_eliminating;i,projected_radial_size-j+1];
%             Index_eliminating=[Index_eliminating;height_size-i+1,projected_radial_size-j+1];
%         end
%     end
% end
% if ~isempty(Index_eliminating)
%     
%     Num_eliminating=size(Index_eliminating,1);
%     Eliminating_list=zeros(1,deviation_angle_size*Num_eliminating);
%     for i=1:Num_eliminating
%         
%         Eliminating_list(deviation_angle_size*(i-1)+1:deviation_angle_size*i)=...
%             ((Index_eliminating(i,1)-1)*projected_radial_size*deviation_angle_size+(Index_eliminating(1,2)-1)*deviation_angle_size+1):...
%             ((Index_eliminating(i,1)-1)*projected_radial_size*deviation_angle_size+(Index_eliminating(1,2))*deviation_angle_size);
%     end
% end
% Reserve_indices=setdiff(1:hist_size,Eliminating_list);
% Hist_out=Hist(:,Reserve_indices);
 end