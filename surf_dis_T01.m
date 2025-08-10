clear; close all;
%A. Meneses-Gutierrez & T. Saito (2025)
%--------------------------------------------------------------------------
% This program estimates surface displacements at GNSS stations
% due to 1 m of slip on triangular subfaults, oriented consistently with
% the uniform coseismic slip direction.

% This program uses the function **TDdispHS** from the software
% developed by Nikkhoo and Walter (2015) to compute displacements (and/or
% strains and stresses) due to triangular dislocations in a half-space.
%
% Reference:
%Nikkhoo, M., & Walter, T. R. (2015). Triangular dislocation: an analytical,
% artefact-free solution. Geophysical Journal International, 201(2), 1117–1139.
%  https://doi.org/10.1093/gji/ggv035]

% This code was developed to produce results presented in:
% "Linking Coseismic Slip and Afterslip in Intraplate Earthquakes:
% A Case Study of the 2016 Central Tottori Earthquake, Japan"
% (2025JB031677).

%Horizontal and vertical component: Z is positive downwards in the code 

filev=append('read_2016Tottori');%observations

fileg=append('faultP_par01');

load(fileg,'trir1','trir2','trir3','dimT');
load(fileg,'disx','disy','disz','dimk');%North, East, UD
nfault=1;
%Slip direction
load('slip_vector_T01','rake','raker');%Slip direction

load(filev,'y_gnss','x_gnss');
load(filev,'x_ref','y_ref','z_ref');
%stations location (y=E; x=N)
obs_y=y_gnss(1:end,1);
obs_x=x_gnss(1:end,1);
obs_z=zeros(length(obs_x),1);
nz0=length(obs_x)-2;%stations installed after the quake (no vertical estimates)

%elastic constants
rigid=3.0E10;
poi=0.25;

basis_x01=zeros;
basis_y01=zeros;
basis_z01=zeros;
%Displacement strike and dip slip components
   %rake=frake;
   for ifault=1:nfault
       ss1=cosd(raker(ifault));
       ds1=sind(raker(ifault));
       ts1=zeros(length(ss1),1);
       %slip components (x,y,z)
       ux_slip=zeros(length(obs_y),dimT(ifault));
       uy_slip=zeros(length(obs_y),dimT(ifault));
       uz_slip=zeros(length(obs_y),dimT(ifault));
       %Conectivity list for triangular elements
       trir=[trir1(1:dimT(ifault),ifault), trir2(1:dimT(ifault),ifault), trir3(1:dimT(ifault),ifault)];
       %Grid points for triangular sub-faults configuration
       dis=[disy(1:dimk(ifault),ifault), disx(1:dimk(ifault),ifault), disz(1:dimk(ifault),ifault)];
       for i=1:length(trir)
           %Angela: GNSS observations (ss and ds unit vector in the slip directions)
           [ue1,un1,uv1] = ...
               TDdispHS(y_ref,x_ref,z_ref*-1.0,...
               [dis(trir(i,1),1) dis(trir(i,1),2) -1*dis(trir(i,1),3)],...
               [dis(trir(i,2),1) dis(trir(i,2),2) -1*dis(trir(i,2),3)],...
               [dis(trir(i,3),1) dis(trir(i,3),2) -1*dis(trir(i,3),3)],...
               ss1,ds1,0,poi);
           ux_ref(:,i)=un1; %% m North
           uy_ref(:,i)=ue1;%east
           uz_ref(:,i)=-uv1;%down
       end %% for each mesh on plate
       %observation points
       for i=1:length(trir)
           %Angela: GNSS observations (ss and ds unit vector in the slip directions)
           [ue1,un1,uv1] = ...
               TDdispHS(obs_y,obs_x,obs_z*-1.0,...
               [dis(trir(i,1),1) dis(trir(i,1),2) -1*dis(trir(i,1),3)],...
               [dis(trir(i,2),1) dis(trir(i,2),2) -1*dis(trir(i,2),3)],...
               [dis(trir(i,3),1) dis(trir(i,3),2) -1*dis(trir(i,3),3)],...
               ss1,ds1,0,poi);
           ux_slip(:,i)=un1-ux_ref(:,i); %% m north
           uy_slip(:,i)=ue1-uy_ref(:,i);%east
           uz_slip(:,i)=-uv1-uz_ref(:,i);%down
       end %% for each mesh on plate
       uz_slip=uz_slip(1:nz0,:);
   end

save surface_disp_T01