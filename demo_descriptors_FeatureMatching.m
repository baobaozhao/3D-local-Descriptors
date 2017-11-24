%  Author: Bao Zhao {zhaobao1988@sjtu.edu.cn}
% this code present the process of generating descriptors on key
% points, and generating PRC curve.


[Model,Model_faces]=read_ply('Data\dragon_vrip_res2.ply');
[Scene,Scene_faces]=read_ply('Data\Scene1_0.5.ply');
NUM_KeyPoints=1000;
% randomly selecting 1000 key points on scene and model respectively.
L_scene=size(Scene,1);
Scene_key_List=zeros(NUM_KeyPoints,1);
for i=1:NUM_KeyPoints
    index=round(rand*L_scene);
    while ~isempty(find(Scene_key_List(1:i)==index, 1))
        index=round(rand*L_scene);
    end
    Scene_key_List(i)=index;
end

Transformation=load('Data\Scene1-rs0.xf');% the true tranformation from model to scene.
t=Transformation(1:3,end);  
R=Transformation(1:3,1:3);
T_Model=(R*Model')';
T_Model=[T_Model(:,1)+t(1),T_Model(:,2)+t(2),T_Model(:,3)+t(3)];
[Model_key_List,~]=knnsearch(T_Model,Scene(Scene_key_List,:),'k',1);


% obtaining mesh resolution.
[~,d1]=knnsearch(Model,Model,'k',2);
mr=mean(d1(:,2));
Support_radius=20*mr;

% calculating LPA.
Model_indices=1:size(Model,1);  Model_LPAs= LPA(Model, Model_indices, 7*mr);
Scene_indices=1:size(Scene,1);  Scene_LPAs= LPA(Scene, Scene_indices, 7*mr);


% calculating LRA/F with the method proposed by GUO et al. "Rotational Projection Statistics for 3D Local Surface Description and Object Recognition"
LRFs_1=GUO_LRF( Model,Model_faces,Model_key_List,mr);
LRFs_2=GUO_LRF( Scene,Scene_faces,Scene_key_List,mr);

Model_LRFs=zeros(NUM_KeyPoints,9);
Scene_LRFs=zeros(NUM_KeyPoints,9);
for i=1:NUM_KeyPoints
   Model_LRFs(i,:)=[LRFs_1{i}(1,:),LRFs_1{i}(2,:),LRFs_1{i}(3,:)]; 
   Scene_LRFs(i,:)=[LRFs_2{i}(1,:),LRFs_2{i}(2,:),LRFs_2{i}(3,:)]; 
end

Model_LRAs=Model_LRFs(:,7:9);
Scene_LRAs=Scene_LRFs(:,7:9);

%%%%%%%%%%%%%%%%%%%%%%%=================== the descriptor of only encoding geometrical information  =====================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model_LD=LD(Model,Model_key_List,Model_LPAs,Support_radius,22,Model_LRAs);
Scene_LD=LD(Scene,Scene_key_List,Scene_LPAs,Support_radius,22,Scene_LRAs);
[Pre_LD,Recall_LD]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_LD,Scene_LD,Support_radius);
LD_Area=generating_AUCpr(Pre_LD,Recall_LD);


%%%%%%%%%%%%%%%%%%%%%%%=================== the descriptors generated with the cartesian coordinate system ==================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model_X_ld=X_ld(Model,Model_key_List,Model_LPAs,Support_radius,5,15,Model_LRFs);
Scene_X_ld=X_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,5,15,Scene_LRFs);
[Pre_X_ld,Recall_X_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_X_ld,Scene_X_ld,Support_radius);
X_ld_Area=generating_AUCpr(Pre_X_ld,Recall_X_ld);

Model_Y_ld=Y_ld(Model,Model_key_List,Model_LPAs,Support_radius,4,14,Model_LRFs);
Scene_Y_ld=Y_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,4,14,Scene_LRFs);
[Pre_Y_ld,Recall_Y_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_Y_ld,Scene_Y_ld,Support_radius);
Y_ld_Area=generating_AUCpr(Pre_Y_ld,Recall_Y_ld);

Model_Z_ld=Z_ld(Model,Model_key_List,Model_LPAs,Support_radius,13,14,Model_LRFs);
Scene_Z_ld=Z_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,13,14,Scene_LRFs);
[Pre_Z_ld,Recall_Z_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_Z_ld,Scene_Z_ld,Support_radius);
Z_ld_Area=generating_AUCpr(Pre_Z_ld,Recall_Z_ld);

Model_Pr_ld=Pr_ld(Model,Model_key_List,Model_LPAs,Support_radius,7,15,Model_LRAs);
Scene_Pr_ld=Pr_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,7,15,Scene_LRAs);
[Pre_Pr_ld,Recall_Pr_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_Pr_ld,Scene_Pr_ld,Support_radius);
Pr_ld_Area=generating_AUCpr(Pre_Pr_ld,Recall_Pr_ld);


Model_XY=XY_spatial(Model,Model_key_List,Support_radius,9,11,Model_LRFs);
Scene_XY=XY_spatial(Scene,Scene_key_List,Support_radius,9,11,Scene_LRFs);
[Pre_XY,Recall_XY]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_XY,Scene_XY,Support_radius);
XY_Area=generating_AUCpr(Pre_XY,Recall_XY);

Model_XZ=XZ_spatial(Model,Model_key_List,Support_radius,6,15,Model_LRFs);
Scene_XZ=XZ_spatial(Scene,Scene_key_List,Support_radius,6,15,Scene_LRFs);
[Pre_XZ,Recall_XZ]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_XZ,Scene_XZ,Support_radius);
XZ_Area=generating_AUCpr(Pre_XZ,Recall_XZ);

Model_YZ=YZ_spatial(Model,Model_key_List,Support_radius,6,17,Model_LRFs);
Scene_YZ=YZ_spatial(Scene,Scene_key_List,Support_radius,6,17,Scene_LRFs);
[Pre_YZ,Recall_YZ]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_YZ,Scene_YZ,Support_radius);
YZ_Area=generating_AUCpr(Pre_YZ,Recall_YZ);


Model_Pr_Z=Pr_Z(Model,Model_key_List,Support_radius,6,15,Model_LRAs);
Scene_Pr_Z=Pr_Z(Scene,Scene_key_List,Support_radius,6,15,Scene_LRAs);
[Pre_Pr_Z,Recall_Pr_Z]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_Pr_Z,Scene_Pr_Z,Support_radius);
Pr_Z_Area=generating_AUCpr(Pre_Pr_Z,Recall_Pr_Z);


Model_XYZ=XYZ_spatial(Model,Model_key_List,Support_radius,8,5,12,Model_LRFs);
Scene_XYZ=XYZ_spatial(Scene,Scene_key_List,Support_radius,8,5,12,Scene_LRFs);
[Pre_XYZ,Recall_XYZ]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_XYZ,Scene_XYZ,Support_radius);
XYZ_Area=generating_AUCpr(Pre_XYZ,Recall_XYZ);


Model_SDASS=SDASS(Model,Model_key_List,Model_LPAs,Support_radius,6,5,14,Model_LRAs);
Scene_SDASS=SDASS(Scene,Scene_key_List,Scene_LPAs,Support_radius,6,5,14,Scene_LRAs);
[Pre_SDASS,Recall_SDASS]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_SDASS,Scene_SDASS,Support_radius);
SDASS_Area=generating_AUCpr(Pre_SDASS,Recall_SDASS);


Model_XY_ld=XY_ld(Model,Model_key_List,Model_LPAs,Support_radius,4,3,12,Model_LRFs);
Scene_XY_ld=XY_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,4,3,12,Scene_LRFs);
[Pre_XY_ld,Recall_XY_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_XY_ld,Scene_XY_ld,Support_radius);
XY_ld_Area=generating_AUCpr(Pre_XY_ld,Recall_XY_ld);



Model_XZ_ld=XZ_ld(Model,Model_key_List,Model_LPAs,Support_radius,4,8,15,Model_LRFs);
Scene_XZ_ld=XZ_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,4,8,15,Scene_LRFs);
[Pre_XZ_ld,Recall_XZ_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_XZ_ld,Scene_XZ_ld,Support_radius);
XZ_ld_Area=generating_AUCpr(Pre_XZ_ld,Recall_XZ_ld);


Model_YZ_ld=YZ_ld(Model,Model_key_List,Model_LPAs,Support_radius,3,8,15,Model_LRFs);
Scene_YZ_ld=YZ_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,3,8,15,Scene_LRFs);
[Pre_YZ_ld,Recall_YZ_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_YZ_ld,Scene_YZ_ld,Support_radius);
YZ_ld_Area=generating_AUCpr(Pre_YZ_ld,Recall_YZ_ld);


Model_XYZ_ld=XYZ_ld(Model,Model_key_List,Model_LPAs,Support_radius,3,3,8,11,Model_LRFs);
Scene_XYZ_ld=XYZ_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,3,3,8,11,Scene_LRFs);
[Pre_XYZ_ld,Recall_XYZ_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_XYZ_ld,Scene_XYZ_ld,Support_radius);
XYZ_ld_Area=generating_AUCpr(Pre_XYZ_ld,Recall_XYZ_ld);


%%%%%%%%%%%%%%%%%%%%%%%=================== the descriptors generated with the polar coordinate system =====================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model_radial_ld=radial_ld(Model,Model_key_List,Model_LPAs,Support_radius,8,14,Model_LRAs);
Scene_radial_ld=radial_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,8,14,Scene_LRAs);
[Pre_radial_ld,Recall_radial_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_radial_ld,Scene_radial_ld,Support_radius);
radial_ld_Area=generating_AUCpr(Pre_radial_ld,Recall_radial_ld);

Model_azimuth_ld=azimuth_ld(Model,Model_key_List,Model_LPAs,Support_radius,5,14,Model_LRFs);
Scene_azimuth_ld=azimuth_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,5,14,Scene_LRFs);
[Pre_azimuth_ld,Recall_azimuth_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_azimuth_ld,Scene_azimuth_ld,Support_radius);
azimuth_ld_Area=generating_AUCpr(Pre_azimuth_ld,Recall_azimuth_ld);

Model_elevation_ld=elevation_ld(Model,Model_key_List,Model_LPAs,Support_radius,11,15,Model_LRFs);
Scene_elevation_ld=elevation_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,11,15,Scene_LRFs);
[Pre_elevation_ld,Recall_elevation_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_elevation_ld,Scene_elevation_ld,Support_radius);
elevation_ld_Area=generating_AUCpr(Pre_elevation_ld,Recall_elevation_ld);

Model_radial_azimuth=radial_azimuth(Model,Model_key_List,Support_radius,6,9,Model_LRFs);
Scene_radial_azimuth=radial_azimuth(Scene,Scene_key_List,Support_radius,6,9,Scene_LRFs);
[Pre_radial_azimuth,Recall_radial_azimuth]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_radial_azimuth,Scene_radial_azimuth,Support_radius);
radial_azimuth_Area=generating_AUCpr(Pre_radial_azimuth,Recall_radial_azimuth);

Model_radial_elevation =radial_elevation (Model,Model_key_List,Support_radius,5,12,Model_LRFs);
Scene_radial_elevation =radial_elevation (Scene,Scene_key_List,Support_radius,5,12,Scene_LRFs);
[Pre_radial_elevation,Recall_radial_elevation]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_radial_elevation ,Scene_radial_elevation ,Support_radius);
radial_elevation_Area=generating_AUCpr(Pre_radial_elevation,Recall_radial_elevation);

Model_azimuth_elevation=azimuth_elevation(Model,Model_key_List,Support_radius,5,12,Model_LRFs);
Scene_azimuth_elevation=azimuth_elevation(Scene,Scene_key_List,Support_radius,5,12,Scene_LRFs);
[Pre_azimuth_elevation,Recall_azimuth_elevation]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_azimuth_elevation,Scene_azimuth_elevation,Support_radius);
azimuth_elevation_Area=generating_AUCpr(Pre_azimuth_elevation,Recall_azimuth_elevation);

Model_radial_azimuth_elevation=radial_azimuth_elevation(Model,Model_key_List,Support_radius,6,4,10,Model_LRFs);
Scene_radial_azimuth_elevation=radial_azimuth_elevation(Scene,Scene_key_List,Support_radius,6,4,10,Scene_LRFs);
[Pre_radial_azimuth_elevation,Recall_radial_azimuth_elevation]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_radial_azimuth_elevation,Scene_radial_azimuth_elevation,Support_radius);
radial_azimuth_elevation_Area=generating_AUCpr(Pre_radial_azimuth_elevation,Recall_radial_azimuth_elevation);

Model_radial_azimuth_ld=radial_azimuth_ld(Model,Model_key_List,Model_LPAs,Support_radius,5,2,14,Model_LRFs);
Scene_radial_azimuth_ld=radial_azimuth_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,5,2,14,Scene_LRFs);
[Pre_radial_azimuth_ld,Recall_radial_azimuth_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_radial_azimuth_ld,Scene_radial_azimuth_ld,Support_radius);
radial_azimuth_ld_Area=generating_AUCpr(Pre_radial_azimuth_ld,Recall_radial_azimuth_ld);

Model_radial_elevation_ld=radial_elevation_ld(Model,Model_key_List,Model_LPAs,Support_radius,4,6,16,Model_LRFs);
Scene_radial_elevation_ld=radial_elevation_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,4,6,16,Scene_LRFs);
[Pre_radial_elevation_ld,Recall_radial_elevation_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_radial_elevation_ld,Scene_radial_elevation_ld,Support_radius);
radial_elevation_ld_Area=generating_AUCpr(Pre_radial_elevation_ld,Recall_radial_elevation_ld);

Model_azimuth_elevation_ld=azimuth_elevation_ld(Model,Model_key_List,Model_LPAs,Support_radius,2,8,13,Model_LRFs);
Scene_azimuth_elevation_ld=azimuth_elevation_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,2,8,13,Scene_LRFs);
[Pre_azimuth_elevation_ld,Recall_azimuth_elevation_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_azimuth_elevation_ld,Scene_azimuth_elevation_ld,Support_radius);
azimuth_elevation_ld_Area=generating_AUCpr(Pre_azimuth_elevation_ld,Recall_azimuth_elevation_ld);

Model_radial_azimuth_elevation_ld=radial_azimuth_elevation_ld(Model,Model_key_List,Model_LPAs,Support_radius,3,2,7,11,Model_LRFs);
Scene_radial_azimuth_elevation_ld=radial_azimuth_elevation_ld(Scene,Scene_key_List,Scene_LPAs,Support_radius,3,2,7,11,Scene_LRFs);
[Pre_radial_azimuth_elevation_ld,Recall_radial_azimuth_elevation_ld]=generate_RPC(T_Model(Model_key_List,:),Scene(Scene_key_List,:),Model_radial_azimuth_elevation_ld,Scene_radial_azimuth_elevation_ld,Support_radius);
radial_azimuth_elevation_ld_Area=generating_AUCpr(Pre_radial_azimuth_elevation_ld,Recall_radial_azimuth_elevation_ld);


%%%%%%%%%%%%%%%%%%%%%%%=================== plot PRC curve =====================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ALL_Area=[X_ld_Area,Y_ld_Area,Z_ld_Area,Pr_ld_Area,XY_Area,XZ_Area,YZ_Area,Pr_Z_Area,...
    XY_ld_Area,XZ_ld_Area,YZ_ld_Area,SDASS_Area,XYZ_Area,XYZ_ld_Area];

[~,List_sort]=sort(ALL_Area,'descend');
Pre_all=[Pre_X_ld';Pre_Y_ld';Pre_Z_ld';Pre_Pr_ld';Pre_XY';Pre_XZ';Pre_YZ';Pre_Pr_Z';...
    Pre_XY_ld';Pre_XZ_ld';Pre_YZ_ld';Pre_SDASS';Pre_XYZ';Pre_XYZ_ld'];
Recall_all=[Recall_X_ld';Recall_Y_ld';Recall_Z_ld';Recall_Pr_ld';Recall_XY';Recall_XZ';Recall_YZ';Recall_Pr_Z';...
    Recall_XY_ld';Recall_XZ_ld';Recall_YZ_ld';Recall_SDASS';Recall_XYZ';Recall_XYZ_ld'];

STR_XYZ={'f_x(ld)','f_y(ld)','f_z(ld)','f_p_r(ld)','f_x_,_y','f_x_,_z','f_y_,_z','SI','f_x_,_y(ld)','f_x_,_z(ld)','f_y_,_z(ld)'...
 ,'SDASS','f_x_,_y_,_z','f_x_,_y_,_z(ld)'};
Mark_style={'-^','-v','-<','->','-p',':o','-s','-+','-d','-x','-h','-o',':^','-*'};
Color=[1,0,0;0,1,0;0,0,1;0,1,1;234/255,122/255,178/255;0,0,0;131/255,67/255,167/255;...
    72/255,106/255,162/255;77/255,157/255,138/255;114/255,183/255,51/255;159/255,117/255,75/255;...
    1,0,1;212/255,150/255,50/255;54/255,142/255,88/255];

figure
hold on
for i=1:14
    STR_XYZ_sort{i}=[STR_XYZ{List_sort(i)},' ','[',num2str(ALL_Area(List_sort(i)),2),']'];
    plot(Pre_all(List_sort(i),:),Recall_all(List_sort(i),:),Mark_style{List_sort(i)},'color',Color(List_sort(i),:),'markersize',5)% 0.2
end

ylim([0,1])
Y=get(gca,'YTick');
str=cell(length(Y),1);
for i=1:1:length(Y)
    str{i}=Y(i);  
end
set(gca,'YTickLabel',str,'YTickMode','manual')

xlim([0,1])
X=get(gca,'XTick');
str=cell(length(X),1);
for i=1:1:length(X)
    str{i}=X(i);   
end
set(gca,'XTickLabel',str,'XTickMode','manual')

ll=legend(STR_XYZ_sort);
set(ll,'Box', 'off','FontSize',8,'Location','northeastoutside')

ylabel('Recall')
xlabel('1-Precision')
grid on
set(gcf,'color',[1,1,1])







ALL_RAE_Area=[LD_Area,radial_ld_Area,azimuth_ld_Area,elevation_ld_Area,radial_azimuth_Area,radial_elevation_Area,azimuth_elevation_Area,radial_azimuth_ld_Area,radial_elevation_ld_Area,...
    azimuth_elevation_ld_Area,radial_azimuth_elevation_Area,radial_azimuth_elevation_ld_Area];

[Sort_zhi,List_sort_RAE]=sort(ALL_RAE_Area,'descend');

STR_RAE={'f(ld)','f_r(ld)','f_a(ld)','f_e(ld)','f_r_,_a','f_r_,_e','f_a_,_e','f_r_,_a(ld)','f_r_,_e(ld)','f_a_,_e(ld)','USC','f_r_,_a_,_e(ld)'};
STR_RAE_sort=cell(12,1);
Pre_all=[Pre_LD';Pre_radial_ld';Pre_azimuth_ld';Pre_elevation_ld';Pre_radial_azimuth';Pre_radial_elevation';Pre_azimuth_elevation';Pre_radial_azimuth_ld';...
   Pre_radial_elevation_ld';Pre_azimuth_elevation_ld';Pre_radial_azimuth_elevation';Pre_radial_azimuth_elevation_ld'];
Recall_all=[Recall_LD';Recall_radial_ld';Recall_azimuth_ld';Recall_elevation_ld';Recall_radial_azimuth';Recall_radial_elevation';Recall_azimuth_elevation';Recall_radial_azimuth_ld';...
   Recall_radial_elevation_ld';Recall_azimuth_elevation_ld';Recall_radial_azimuth_elevation';Recall_radial_azimuth_elevation_ld'];
figure
hold on
for i=1:12
    STR_RAE_sort{i}=[STR_RAE{List_sort_RAE(i)},' ','[',num2str(ALL_RAE_Area(List_sort_RAE(i)),2),']'];
    plot(Pre_all(List_sort_RAE(i),:),Recall_all(List_sort_RAE(i),:),Mark_style{List_sort_RAE(i)},'color',Color(List_sort_RAE(i),:),'markersize',5)% 0.2
end

ylim([0,1])
Y=get(gca,'YTick');
str=cell(length(Y),1);
for i=1:1:length(Y)
    str{i}=Y(i);   
end
set(gca,'YTickLabel',str,'YTickMode','manual')

xlim([0,1])
X=get(gca,'XTick');
str=cell(length(X),1);
for i=1:1:length(X)
    str{i}=X(i);   
end
set(gca,'XTickLabel',str,'XTickMode','manual')


set(gca,'FontSize',8)
ll=legend(STR_RAE_sort);
set(ll,'Box', 'off','FontSize',8,'Location','northeastoutside')
ylabel('Recall')
xlabel('1-Precision')
grid on
set(gcf,'color',[1,1,1])

