clc, clear, close all;
% Zhao_4 SHPSG function

% VARIABLES TO MANIPULATE
EI = 1; % elongation index
FI = 0.5; % flatness index
D2_8 = 0.6;
D9_15 = 0;

coeff = SHPSG(EI,FI,D2_8,D9_15);
theta = linspace(0,pi,50);
phi = linspace(0,2*pi,50);
[theta,phi] = meshgrid(theta,phi);

% view particle in matlab
figure;
view(3);
[x,y,z] = sph2car(coeff, phi, theta);
surf(x,y,z,'EdgeColor','none','FaceColor','black','FaceAlpha',0.4);
axis equal;
light

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xyz = [x(:), y(:), z(:)];

% matrix dimensions
dim = 100;
[Xlim,Ylim,Zlim] = generateBox(x,y,z);

ix = 2*linspace(Xlim(1),Xlim(2),dim);
iy = 2*linspace(Ylim(1),Ylim(2),dim);
iz = 2*linspace(Zlim(1),Zlim(2),dim);

% create a matrix with all the valid test points
testpts = zeros(dim*dim*dim,3);
count = 1;
for x=1:dim
    for y=1:dim
        for z=1:dim
            testpts(count,1) = ix(x);
            testpts(count,2) = iy(y);
            testpts(count,3) = iz(z);
            count = count+1;
        end
    end
end

inspherical(:,1)=sqrt(testpts(:,1).^2+testpts(:,2).^2+testpts(:,3).^2); %r


for i=1:length(inspherical)

    %%theta 0 to pi

if testpts(i,3)>0
    inspherical(i,2) = atan(sqrt(testpts(i,1)^2+testpts(i,2)^2)/testpts(i,3));
elseif testpts(i,3)<0
    inspherical(i,2) = pi+atan(sqrt(testpts(i,1)^2+testpts(i,2)^2)/testpts(i,3));
else
 inspherical(i,2)=(pi/2);
end

%% phi 0 to 2pi

if testpts(i,1)>0
inspherical(i,3)=atan(testpts(i,2)/testpts(i,1));
elseif testpts(i,1)<0 && testpts(i,2)>0
inspherical(i,3)=pi+atan(testpts(i,2)/testpts(i,1));
elseif testpts(i,1)<0 && testpts(i,2)<0
inspherical(i,3)=-pi+atan(testpts(i,2)/testpts(i,1));
elseif testpts(i,1)==0 && testpts(i,2)>0
inspherical(i,3)=pi/2;
elseif testpts(i,1)==0 && testpts(i,2)<0
inspherical(i,3)=-pi/2;
end

end

inspherical(:,3)=pi+inspherical(:,3);

[px,py,pz] = sph2car(coeff, inspherical(:,3), inspherical(:,2));

r=sqrt(px.^2+py.^2+pz.^2);


% represent particle using 3D matrix
img = zeros(dim,dim,dim);
count = 1;

% match the indexes of points within the object with the corresponding index on the matrix 
for i=1:length(r)
    if inspherical(i,1)>r(i)
    img(count) = 0;
    else
    img(count) = 1;
    end
    count = count+1;
end

% save particle
outFile = 'particle.tif';
imwrite(img(:,:,1),outFile);
for ii=2:size(img,3)
    imwrite(img(:,:,ii),outFile,'WriteMode','append')
end




function [Xlim,Ylim,Zlim] = generateBox(x,y,z)

    % generate radiograph 
    Xlim=[min(x,[],"all"); max(x,[],"all")];
    Ylim=[min(y,[],"all"); max(y,[],"all")];
    Zlim=[min(z,[],"all"); max(z,[],"all")];
    
    % create edges
    X=[Xlim(1);Xlim(2);Xlim(2);Xlim(1);Xlim(1);Xlim(2);Xlim(2);Xlim(1)];
    Y=[Ylim(1);Ylim(1);Ylim(1);Ylim(1);Ylim(2);Ylim(2);Ylim(2);Ylim(2)];
    Z=[Zlim(1);Zlim(1);Zlim(2);Zlim(2);Zlim(1);Zlim(1);Zlim(2);Zlim(2)];
    
    Xfront=[X(1);X(2);X(3);X(4)];Xback=[X(5);X(6);X(7);X(8)];Xtop=[X(4);X(3);X(7);X(8)];Xbottom=[X(1);X(2);X(6);X(5)];Xleft=[Y(1);Y(5);Y(4);Y(8)];Xright=[Y(2);Y(6);Y(7);Y(3)];
    
    Yfront=[Y(1);Y(2);Y(3);Y(4)];Yback=[Y(5);Y(6);Y(7);Y(8)];Ytop=[Y(4);Y(3);Y(7);Y(8)];Ybottom=[Y(1);Y(2);Y(6);Y(5)];Yleft=[Y(1);Y(5);Y(4);Y(8)];Yright=[Y(2);Y(6);Y(7);Y(3)];
    
    Zfront=[Z(1);Z(2);Z(3);Z(4)];Zback=[Z(5);Z(6);Z(7);Z(8)];Ztop=[Z(4);Z(3);Z(7);Z(8)];Zbottom=[Z(1);Z(2);Z(6);Z(5)];
    
    hold on
    fill3(Xfront,Yfront,Zfront,'r','FaceColor','none','EdgeColor','r','LineWidth',3)
    fill3(Xback,Yback,Zback,'r','FaceColor','none','EdgeColor','r','LineWidth',3)
    
    fill3(Xtop,Ytop,Ztop,'r','FaceColor','none','EdgeColor','r','LineWidth',3)
    fill3(Xbottom,Ybottom,Zbottom,'r','FaceColor','none','EdgeColor','r','LineWidth',3)
end

function [fvec] = SHPSG(EI,FI,D2_8,D9_15)

    Fvec = zeros(4,3);

    % determine C0 and C1 with a unit max principal dimension
    arr = [0 0 0;
        -1 1j 0;
        0 0 sqrt(2);
        1 1j 0];

    fvec_sphere = -sqrt(pi/6)*arr; % a sphere with a unit diameter
    Fvec(:,1) = fvec_sphere(:,1);
    Fvec(:,2) = EI*fvec_sphere(:,2);
    Fvec(:,3) = (FI*EI)*fvec_sphere(:,3);

    d1 = norm(Fvec,'fro');

    % Determine D2 and D9
    alpha = 1.387;
    beta = 1.426;

    D_2 = D2_8/(1+(2/3)^alpha+(2/4)^alpha+(2/5)^alpha+(2/6)^alpha+(2/7)^alpha+(2/8)^alpha)*d1;
    D_9 = D9_15/(1+(9/10)^beta+(9/11)^beta+(9/12)^beta+(9/13)^beta+(9/14)^beta+(9/15)^beta)*d1;

    % Determine D3-D8 and D10-D15
    % assume all descriptors have three identical decomposition at x, y, and z-axis
    I = zeros(16,3);
    I(2,:) = [1 1 1];
    I(3,:) = [D_2 D_2 D_2];

    for c = 3:8
        dn = D_2*((c-1)/2)^(-alpha)/sqrt(3);
        I(c+1,:) = [dn dn dn];
    end

    for c = 9:14
        dn = D_9*((c-1)/9)^(-beta)/sqrt(3); % or should it be c-1
        I(c+1,:) = [dn dn dn];
    end

    L = zeros(16^2, 3);
    N = zeros(16^2, 3);


    % Randomly generate P
    % Randomly generate P including C1'-C15' with c_n^(-m)=(-1)^m*c_n^m*

    for n = 1:15
        J = ones(n+1,3)-2*rand(n+1,3); % [-1,1]
        K = flipud(J(1:n,:));
        i = (-1)^n; A = [i i i];
        B = K.*A;
        L(n^2+1:(n+1)^2,:) = [J; B];
        M = ones(n,3)-2*rand(n,3);
        M1 = [M; 0 0 0];
        N(n^2+1:(n+1)^2,:) = [M1; flipud((-1)^(n+1)*M)];
    end

    P = L+N*1j;
    P(1,:) = ones(1,3)-2*rand(1,3);

    % Calculate D1-D16 with the SH coefficients of P
    Q = conj(P);
    R = zeros(16,3);

    for n = 0:15
        a = n^2+1;
        b = (n+1)^2; % matlab end is inclusive but python is not
        R(n+1,:) = sqrt(real(sum(P(a:b,:).*Q(a:b,:))));
    end

    % Determine C2-C15 by making the descriptors of P equal to D2-D15
    fvec = zeros(16^2,3);
    fvec(1:4,:) = Fvec;

    for d = 3:16
        a = (d-1)^2+1;
        b = d^2;

        for n = 1:3
            fvec(a:b,n) = P(a:b,n)/R(d,n)*I(d,n);
        end
    end
end

% convert spherical to cartesian coordinates with SH expansion
function [x,y,z] = sph2car(coeff, phi, theta)
    
    x = 0; y = 0; z = 0;
    index = 1;

    % compute spherical harmonic of degree N and order M
    for n = 0:8
        order = -n:n;

        for m = 0:2*n
            Y = harmonicY(n,order(m+1),theta,phi);
            x = x + coeff(index,1)*Y;
            y = y + coeff(index,2)*Y;
            z = z + coeff(index,3)*Y;
            index = index + 1;
        end
    end

    x = real(x);
    y = real(y);
    z = real(z);
end
