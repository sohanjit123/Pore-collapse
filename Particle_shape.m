clc,clear all,close all
% Image directory (change this to segment others)
imageDirectory = '/Users/sohanjitghosh/Desktop/OPTIMISATION/';

imageName = '77.tif';


%% Load tif image
imagePath =[imageDirectory,imageName];

Itmp = imread(imagePath); % Read first frame to get partial size
info = imfinfo(imagePath); % Read information to get length
numstack = length(info); % Get length / number of slices
Img = zeros(size(Itmp,1),size(Itmp,2),numstack); % Construct empty image
for ii = 1:numstack % Read image by slice
    currentImage = imread(imagePath,ii,'Info',info);
    Img(:,:,ii) = currentImage;
end
% disp(['Finished loading image ',imageName]);


% %% Normal Marching Cubes
% 
% isovalue=0.05;
% 
% [faces,verts] = extractIsosurface(Img,isovalue);
% % [faces,verts] = isosurface(Img,isovalue);
% 
% 
% figure
% p = patch(Faces=faces,Vertices=verts);
% isonormals(Img,p)
% view(3)
% set(p,FaceColor=[0.5 0.5 0.5])
% set(p,EdgeColor="none")
% camlight
% lighting flat

%% Generalised Marching Cubes with Gaussian Filter

Img(Img>0)=1;

Img=smooth3(Img,'gaussian',5);

s = regionprops3(Img,"all");

isovalue=0.05;

[faces,verts] = extractIsosurface(Img,isovalue);
% [faces,verts] = isosurface(Img,isovalue);

% figure(1)
 p = patch(Faces=faces,Vertices=verts);
isonormals(Img,p)
view(3)
set(p,FaceColor=[0.5 0.5 0.5])
% set(p,FaceAlpha=0.1)
set(p,EdgeColor="none")
camlight
lighting flat

x=double(verts(:,1));
y=double(verts(:,2));
z=double(verts(:,3));


% [k,av] = convhull(x,y,z);


C_x=s.Volume/s.ConvexVolume;

S=((36*pi*s.Volume*s.Volume)^(1/3))/s.SurfaceArea;

PA=sort(s.PrincipalAxisLength,"descend");
a=PA(1);
b=PA(2);
c=PA(3);


EI=b/a;
FI=c/b;


N=1500;

figure(2)
S_v=reducepatch(p,N+1);

P1=patch('Faces',S_v.faces,'Vertices',S_v.vertices);

view(3)
set(P1,FaceColor=[0.5 0.5 0.5])
set(P1,EdgeColor="none")
set(P1,FaceAlpha=0.1)
camlight
lighting flat

%%Curvature of the vertices

  [K_m,K_G,Dir1,Dir2,K1,K2]=patchcurvature(S_v,true);

  % figure, title('Principal A');
  %   p1=S_v.vertices-2*Dir1; p2=S_v.vertices+2*Dir1;       
  %   plot3([p1(:,1) p2(:,1)]',[p1(:,2) p2(:,2)]',[p1(:,3) p2(:,3)]','g-');
  %   axis equal; view(3) 

  D=bwdist(~Img);
  D=double(D);

  rad=  max(D,[],"all");

for i=1:size(D,1)
    for j=1:size(D,2)
        for k=1:size(D,3)

            if D(i,j,k)==rad
                center=[i j k];
                break;
            end

        end
    end
end

% [V,F]=icosphere(4);
% 
% V=V*rad+center.*ones(size(V));
% pins = patch(Faces=F,Vertices=V);
% 
% S_inscir=reducepatch(pins,N+1);

%% Get projection onto maximum inscribed sphere

Numvertices=length(S_v.vertices);

hold on

for i=1:Numvertices

vec=S_v.vertices(i,:)-center;
l=vec(1)/norm(vec);
m=vec(2)/norm(vec);
n=vec(3)/norm(vec);


new_vec(i,:)=rad*[l m n]+center;
% scatter3(new_vec(i,1),new_vec(i,2),new_vec(i,3),'r.')
% hold on


end

S_ins=S_v;
S_ins.vertices=new_vec;

P2=patch('Faces',S_ins.faces,'Vertices',S_ins.vertices);
set(P2,FaceAlpha=0.5)

[F_m,F_G,Fir1,Fir2,F1,F2]=patchcurvature(S_ins,true);

%%  particle Roundess calculation

Numelements=length(S_v.faces);

num=0;den=0;

for i=1:1

FP=S_v.vertices(S_v.faces(i,1),:);
SP=S_v.vertices(S_v.faces(i,2),:);
TP=S_v.vertices(S_v.faces(i,3),:);


L1=norm(FP-SP);
L2=norm(SP-TP);
L3=norm(FP-TP);

s = (L1+L2+L3)/2;
area = sqrt(s*(s-L1)*(s-L2)*(s-L3));


LR=((F_m(S_v.faces(i,1))/K_m(S_v.faces(i,1)))+(F_m(S_v.faces(i,2))/K_m(S_v.faces(i,2)))...
+(F_m(S_v.faces(i,3))/K_m(S_v.faces(i,3))))/3;


num=num+area*LR;

den=den+area;
end

R_M=num/den;
