clear; close all;
%A. Meneses-Gutierrez & T. Saito (2025)
%--------------------------------------------------------------------------
% This program estimates the normal vector, strike direction, and dip
% direction for each triangular subfault, in (x, y, z) coordinates aligned
% with (North, East, down). The slip vector is defined relative to the strike
% direction of a corresponding rectangular fault.
%
% This code was developed to produce results presented in:
% "Linking Coseismic Slip and Afterslip in Intraplate Earthquakes:
% A Case Study of the 2016 Central Tottori Earthquake, Japan"
% (2025JB031677).

%dislocation by the slip and rake
%  ss : strike slip 
%  ds : dip slip  
%  ts : tensile slip
%output vx,vy,vz;n_vx,n_vy,n_vz;ed;es;n
%This program consider one fault plane
%Uniform slip input
load('faultP_par01',...
    'trir1','trir2','trir3','dimT','X','Y','Z',...
    'disy','disx','disz','dimk','ICrx','ICry');
%Rake: slip direction
% Assumed slip direction for calculation
raker=-2;
%Delaunay triangulation: trir1 trir2 trir3 (dimT)
%Centroids for each element of the triangular patches: X, Y, Z
%Dislocation vertex: disx, disy, disz (dimk)
%Projection in local coordinate system (along strike and along the dip): ICrx, ICry
%Reshape data
Xr=zeros;
Yr=zeros;
Zr=zeros;
slipr=zeros;
y_cal=zeros;
z_cal=zeros;
x_cal=zeros;
x_ref=zeros;
y_ref=zeros;
z_ref=zeros;
s_ref=zeros;

%Dip and strike slip component (uniform case)
%1 m slip
ss0=zeros(length(X),1)+cosd(raker);
ds0=zeros(length(X),1)+sind(raker);
ts0=zeros(length(X),1);
rake=zeros(length(X),1)+raker;

%estimate slip direction vectors and normal vectors
vx=zeros;
vy=zeros;
vz=zeros;
ed=zeros;
es=zeros;
disr=zeros;
% Subfault configuration
nf=0;
trir=[trir1(1:dimT,1), trir2(1:dimT,1),trir3(1:dimT,1)];
dis=[disy(1:dimk,1), disx(1:dimk,1), disz(1:dimk,1)];
%
n=zeros;
n_vx=zeros;
n_vy=zeros;
n_vz=zeros;
for i=1:length(trir(:,1))
    nf=nf+1;
    y(1)=dis((trir(i,1)),1);
    x(1)=dis((trir(i,1)),2);
    z(1)=dis((trir(i,1)),3);
    y(2)=dis((trir(i,2)),1);
    x(2)=dis((trir(i,2)),2);
    z(2)=dis((trir(i,2)),3);
    y(3)=dis((trir(i,3)),1);
    x(3)=dis((trir(i,3)),2);
    z(3)=dis((trir(i,3)),3);
    %
    normVec = cross([x(2);y(2);z(2)]-[x(1);y(1);z(1)], [x(3);y(3);z(3)]-[x(1);y(1);z(1)]);
    normVec = normVec./norm(normVec);
    %
    strikeV                    = -[-sin(atan2(normVec(2),normVec(1))) cos(atan2(normVec(2),normVec(1))) 0];
    dipV                       = cross(normVec,strikeV);
    slipC                     = [ss0(nf,1) ds0(nf,1) ts0(nf,1)];
    slipVec                      = [strikeV(:) dipV(:) normVec(:)] * slipC(:);
    %Slip direction vector (x,y,z)
    vx(nf,1)=slipVec(1,1);
    vy(nf,1)=slipVec(2,1);
    vz(nf,1)=slipVec(3,1);
    %
    ed(nf,1)=dipV(1,1);
    ed(nf,2)=dipV(1,2);
    ed(nf,3)=dipV(1,3);
    %
    es(nf,1)=strikeV(1,1);
    es(nf,2)=strikeV(1,2);
    es(nf,3)=strikeV(1,3);
    % vector normal to each sub-fault surface
    n(nf,1)=normVec(1,1);
    n(nf,2)=normVec(2,1);
    n(nf,3)=normVec(3,1);
    %
end
%
%Visual test for slip vectors
figure('Name','normal vector,n','Position',[1500 0 500 350]);
trisurf(trir,dis(:,1),dis(:,2),dis(:,3),'FaceColor','c');hold on;
%end
plot3(Y,X,Z,'*r');
quiver3(Y,X,Z,n(:,2),n(:,1),n(:,3));
xlabel('Y, East [km]');ylabel('X, North [km]');
set(gca,'ZDir','reverse');

set(gca,'ZDir','reverse');
xlabel('Y, East [km]');ylabel('X, North [km]');

save('slip_vector_T01', '-v7.3')

caltime=toc

