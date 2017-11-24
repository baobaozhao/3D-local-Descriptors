function LocalFrames=GUO_LRF( Model,Model_faces,Model_List,mr)
Model_mesh.vertices=Model;
Model_mesh.faces=Model_faces;
L_Model_faces=length(Model_faces);
Model_areas=zeros(L_Model_faces,1);
for i=1:length(Model_faces)
   
    V1=Model(Model_faces(i,1),:)-Model(Model_faces(i,2),:);
    V2=Model(Model_faces(i,1),:)-Model(Model_faces(i,3),:);
    Model_areas(i)=norm(cross(V1,V2))/2;
    
end
Model_mesh.FaceArea=Model_areas;
Model_mesh.FaceCenter=[(Model(Model_faces(:,1),1)+Model(Model_faces(:,2),1)+Model(Model_faces(:,3),1))/3,...
    (Model(Model_faces(:,1),2)+Model(Model_faces(:,2),2)+Model(Model_faces(:,3),2))/3,(Model(Model_faces(:,1),3)+Model(Model_faces(:,2),3)+Model(Model_faces(:,3),3))/3];
Model_mesh.keypntIdx=Model_List;
% LocalFrames = LRFforMeshFuncT(Model_mesh, Model_List, mr,mr);
% TriSI = TriSIFunc(Model_mesh,mr, mr/20, 15, LocalFrames);
 LocalFrames = LRFforMeshFuncT(Model_mesh, Model_List, 20*mr,5*mr);
end

