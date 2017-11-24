function LFs = LRFforMeshFuncT(mesh, keypntIdx, neighborSize,edgeThreshold)
% LRFs = LRFforMeshFunc(mesh, keypntIdx, neighborSize, edgeThreshold)

% Author : Yulan Guo {yulan.guo@nudt.edu.cn}
%  College of Electronics Science and Engineering, National University of
%  Defense Technology
%
% This function takes a mesh as an input, and generate local reference frames (LRFs) for a set of
% keypoints as an output.
%
% Arguments : mesh - with vertices and faces           
%                       keypntIdx - the indices of keypoints on a mesh
%                       neighborSize - the size of neighborhood to define a
%                                                  local surface for a selected keypoint
%                       edgeThreshold - the threshold which is used to
%                       exclude the faces with too long edges
% Return :         LRFs - local reference frames (LRFs) corresponding to all
%                                     keypoints
% Copyright : This code is written by Yulan Guo {yulan.guo@nudt.edu.cn},
%              College of Electronics Science and Engineering, National University of
%              Defense Technology. The code CANNOT be distributed to any
%              third party without the explicit permission of Yulan Guo.
%              Under permission, the code may be used  for research purposes with
%              citation to the following references.
%References:
%           [1] Y Guo, F Sohel, M Bennamoun, M Lu, J Wan. TriSI: A distinctive local surface descriptor for 3D modeling and object recognition
%               8th International Conference on Computer Graphics Theory and Applications (GRAPP),2013
%           [2] Y Guo, F Sohel, M Bennamoun, J Wan, M Lu. A novel local surface feature for 3D object recognition under clutter and occlusion.
%               Information Sciences 293, 196-213, 2015
%
% Disclaimer : This code is provided as is without any warranty.

BucketSize = floor(length(mesh.vertices)/100);
kdtreeFaceCenter = KDTreeSearcher(mesh.FaceCenter,'Distance','euclidean','BucketSize',BucketSize);

for i = 1:length(keypntIdx)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    keypnt = mesh.vertices(keypntIdx(i),:);
    [neighborFaceIdx,neighborFaceDis]= rangesearch(kdtreeFaceCenter,keypnt,neighborSize);    
    neighborFaceIdx = cell2mat(neighborFaceIdx);
    neighborFaceDis = cell2mat(neighborFaceDis);
    neighborFaceLength = length(neighborFaceIdx);
    Len1 = floor(neighborFaceLength*0.8);
    neighborIdx = [];    
    M = zeros(3,3);
    totalArea = 0;
    tiotalDisW = 0;
    faceValid = ones(neighborFaceLength,1);
    for j=1:neighborFaceLength
        vertIdx = mesh.faces(neighborFaceIdx(j),:);  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculate the length of each edge of the face
        edges(1) = norm(mesh.vertices(vertIdx(1),:)-mesh.vertices(vertIdx(2),:));
        edges(2) = norm(mesh.vertices(vertIdx(3),:)-mesh.vertices(vertIdx(3),:));
        edges(3) = norm(mesh.vertices(vertIdx(3),:)-mesh.vertices(vertIdx(1),:));    
        maxEdge = max(edges);    
        if maxEdge>edgeThreshold
            faceValid(j) = 0;
            continue;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if j<Len1%Inner neighboring faces
            neighborIdx = [neighborIdx,vertIdx];
            area(j) = mesh.FaceArea(neighborFaceIdx(j));
            totalArea = totalArea+area(j);            
            disW(j) = (neighborSize-neighborFaceDis(j))^2;
            tiotalDisW = tiotalDisW+disW(j);
            center_i = zeros(3,3);
            for ii=1:3
                temp11 = mesh.vertices(vertIdx(ii),:)-keypnt;
                for jj=1:3
                    center_i = center_i+temp11'*(mesh.vertices(vertIdx(jj),:)-keypnt);
                end
                center_i = center_i+temp11'*temp11;
                displace(ii,:) = temp11;
            end
            displaceTotal{j} = displace;
            M = M+center_i*area(j)*disW(j);
        else%outer neighboring faces
            dis1(j) = norm(keypnt-mesh.vertices(vertIdx(1),:));
            dis2(j) = norm(keypnt-mesh.vertices(vertIdx(2),:));
            dis3(j) = norm(keypnt-mesh.vertices(vertIdx(3),:));    
            if dis1(j)<neighborSize && dis2(j)<neighborSize && dis3(j)<neighborSize
                neighborIdx = [neighborIdx,vertIdx];
                area(j) = mesh.FaceArea(neighborFaceIdx(j));
                totalArea = totalArea+area(j);            
                disW(j) = (neighborSize-neighborFaceDis(j))^2;
                tiotalDisW = tiotalDisW+disW(j);
                center_i = zeros(3,3);
                for ii=1:3
                    temp11 = mesh.vertices(vertIdx(ii),:)-keypnt;
                    for jj=1:3
                        center_i = center_i+temp11'*(mesh.vertices(vertIdx(jj),:)-keypnt);
                    end
                    center_i = center_i+temp11'*temp11;
                    displace(ii,:) = temp11;
                end
                displaceTotal{j} = displace;   
                M = M+center_i*area(j)*disW(j);
            end
        end          
    end
    
    M = M/totalArea/tiotalDisW;    
    %%%%%%
    if isnan(M(1,1)) ==1 || isinf(M(1,1)) ==1
        LFs{i,1} = eye(3,3);
        continue;
    end
    %%%%%
    [V,D] = eig(M);
    lamda = [D(1,1),D(2,2),D(3,3)];
    [temp, idxMinLam] = min(lamda);
    [temp, idxMaxLam] = max(lamda);
    xtemp = V(:,idxMaxLam);
    ztemp = V(:,idxMinLam);
    xDisplace = 0;
    zDisplace = 0;

    %disambiguating the sign of x- y- and z- axis
    for j=1:neighborFaceLength
        if faceValid(j) == 0    %invalid face
            continue;
        end
        if j>=Len1
            rf = dis1(j)<neighborSize && dis2(j)<neighborSize && dis3(j)<neighborSize;
            if rf==0
                continue;
            end
        end
        displace = displaceTotal{j};   
        for ii=1:3
            
            xDisplace = xDisplace + displace(ii,:)*xtemp*area(j)*disW(j);
            zDisplace = zDisplace + displace(ii,:)*ztemp*area(j)*disW(j);
        end      
    end
    if xDisplace>0
        xAxis = xtemp;
    else
        xAxis = -xtemp;
    end
   if zDisplace>0
        zAxis = ztemp;
    else
        zAxis = -ztemp;        
    end 
    yAxis = cross(zAxis,xAxis);    
    %get the local coordinates of neighboring points
    rotation = [xAxis';yAxis';zAxis'];
    LFs{i,1} = rotation;
end