      %A. Meneses-Gutierrez & T. Saito (2025)
      %--------------------------------------------------------------------------
      
      % This program performs the stress change inversion using the results from
      % the previous subprograms
      % (Afterslip is given in m)

      % This code was developed to produce results presented in:
      % "Linking Coseismic Slip and Afterslip in Intraplate Earthquakes:
      % A Case Study of the 2016 Central Tottori Earthquake, Japan"
      % (2025JB031677).

    
      clear; close all;
    
      load('faultP_par01','disy','disx','disz','trir1','trir2','trir3',...
          'dimT','dimk','ICrx','ICry');%for figure
      load('faultP_par01','ICr','Rot1','Rot2','Rot3');%dislocations
      load('faultP_par01','fleng','fwidth');
      load('area_pos_B01','x_aft','z_aft');%%basis functions centers
      load('area_pos_B01','basisf','grid_y_aft','grid_x_aft','grid_z_aft');
      load('slip_vector_T01','rake');%slip direction
      %data
      filev=append('read_2016Tottori');
      load(filev,'y_gnss','x_gnss','p_ns','p_ew','p_ud','pe_ew',...
          'pe_ns','pe_ud','lonA','latA');
      obs_vy_r=p_ew(1:end,1);
      obs_vx_r=p_ns(1:end,1);
      obs_vz_r=-p_ud(1:end-2,1);%Z positive downwards
      nz0=length(obs_vz_r);
      err_vy_r=pe_ew(1:end,1);
      err_vx_r=pe_ns(1:end,1);
      err_vz_r=pe_ud(1:end-2,1);
      stlon=lonA(1:end,1);
      stlat=latA(1:end,1);
      % response function
      load('surface_disp_T01','obs_x','obs_y');
      load('stressbase_T01','uy_basis','ux_basis','uz_basis');
      load('stressbase_T01','slipG_f0');

      %elastic constants
      rigid=3.0E10;
      poi=0.25;
      %Damping parameter
      alpha=10.0;
      %Area of triangle: Heron Formula
      pitag=sqrt(fleng.^2+fwidth.^2);
      semip=(pitag+fleng+fwidth)./2.0;
      area1=sqrt(semip.*(semip-fwidth).*(semip-fleng).*(semip-pitag));
    
      %NN constraints (NN=1)
      NN=0;
      %Horizontal and vertical component  horver=1; only horizontal=0
      horver=1;

      opath=append('results');
      mkdir(opath);

      name_mod='after_drop';
      %%
      trir=[trir1(1:dimT,1), trir2(1:dimT,1),trir3(1:dimT,1)];
      dis=[disy(1:dimk,1), disx(1:dimk,1), disz(1:dimk,1)];

      surf_y= obs_vy_r; %m
      surf_x= obs_vx_r;
      surf_z= obs_vz_r;
    
      %--------------------
      % Inversion

      %Horizontal and vertical components consideration
      if(horver==1)
          d1=[surf_x; surf_y; surf_z];
          w1=[1./err_vx_r; 1./err_vy_r; 1./err_vz_r];
          wori=[err_vx_r; err_vy_r; err_vz_r];
          Wt = diag(w1);% W^{1/2}, for weighted least squares
          dw=Wt*d1;
          Gall=[ux_basis; uy_basis; uz_basis];
          GG=Wt*Gall;
      else
          d1=[surf_x; surf_y];
          w1=[1./err_vx_r; 1./err_vy_r];
          wori=[err_vx_r; err_vy_r];
          Wt = diag(w1);
          dw=Wt*d1;
          Gall=[ux_basis; uy_basis];
          GG=Wt*Gall;
      end
      size1=size(uy_basis);
      num_model=size1(2); %number of model parameters 129
      %Damping
      G1=eye(num_model,num_model);

      sa1=zeros;
      rms1=zeros;
      maxslip=zeros;
      minslip=zeros;
      wrms1=zeros;
      semiL=zeros;
    
    
      dinv=[dw; zeros(num_model,1)];
      Ginv=[GG; alpha*G1];

      %Angela
      if(NN==1)
          mm = lsqnonneg(Ginv,dinv);
      else
          mm =lsqminnorm(Ginv,dinv);
      end

      %Weighted solution

      cal_all=Gall*mm;
      ncal=length(cal_all);
      cal_vx=cal_all(1:length(surf_x));
      cal_vy=cal_all(length(surf_x)+1:length(surf_x)+length(surf_y));
      Resi_vx=surf_x-cal_vx;
      Resi_vy=surf_y-cal_vy;
      if(horver==1)
          cal_vz= cal_all(length(surf_x)+length(surf_y)+1:end);
          Resi_vz=surf_z-cal_vz;
      end


      diff01=d1-cal_all;

      %least square solution
      sa1=(dinv-Ginv*mm)'*(dinv-Ginv*mm);

      totalstress=basisf*mm;
      totalslip=slipG_f0*mm;

      maxslip=max(abs(totalslip));
      minslip=min((totalslip));
      %Spatial distribution (draft)
      text1=append('alpha = ',num2str(alpha));
      hfig=figure('doublebuffer','off','Visible','Off');
      subplot(2,2,1);
      s=scatter(ICrx(1:dimT(1),1),ICry(1:dimT(1),1),8,...
          totalstress,'filled'); hold on;title(text1);
      axis equal;axis([-9 21 -30 2]);
      shading flat;
      d=colorbar;d.Label.String = 'stress [MPa]';
      clim([-1.5 1.5]);

      subplot(2,2,2);
      s=scatter(ICrx(1:dimT(1),1),ICry(1:dimT(1),1),8,...
          totalslip,'filled'); hold on;
      s.SizeData = 50;
      %Basis functions for afterslip inversion
      plot(x_aft,z_aft,'xk');
      xlabel('Y, East [km]');ylabel('X, North [km]');zlabel('Z [km]');
      d=colorbar;
      clim([-0.3 0.3]);
      d.Label.String = 'Recovered afterslip [m]';
      axis equal;axis([-9 21 -30 2]);
      savefig=[sprintf('results/Inversion_%s_%s.pdf',name_mod,num2str(alpha))];
      saveas(hfig,savefig);

      %Observation vs calculation for each model
      hfig=figure('doublebuffer','off','Visible','Off');
      if(horver==1)
          subplot(2,1,1);
      end
      scale1=350;
      plot(grid_y_aft,grid_x_aft,'.','Color','k'); hold on;
      quiver(obs_y,obs_x,scale1*surf_y,scale1*surf_x,'off',"k");
      quiver(402,560,0,scale1*0.05,'off',"b");
      text1=append('','5cm Obs.'); text(400,560,text1,'FontSize',9);
      %Model
      quiver(obs_y,obs_x,scale1*cal_vy,scale1*cal_vx,'off',"r");
      axis equal;axis([280 420 550 640]);
      if(horver || 0)
          subplot(2,1,2);
          quiver(obs_y(1:length(surf_z)),obs_x(1:length(surf_z)),zeros(length(surf_z),1),-scale1*surf_z,'off',"b");
          hold on;
          quiver(402,560,0,scale1*0.1,'off',"b");
          text1=append('','10cm Obs.'); text(400,560,text1,'FontSize',9);
          quiver(obs_y(1:length(surf_z)),obs_x(1:length(surf_z)),zeros(length(surf_z),1),-scale1*cal_vz,'off',"r");
          axis equal;axis([280 420 550 640]);
      end
      savefig=[sprintf('results/Displacement%s_%s.pdf',name_mod,num2str(alpha))];
      saveas(hfig,savefig);
    
      %output shear traction
      text1=append('results/Traction_inv_',name_mod,'.sd');
      fid02=fopen(text1,'w');
      for ind=1:dimT
          fprintf(fid02,'%8.4f %8.4f %8.4f \n',ICrx(ind,1),...
              ICry(ind,1),totalstress(ind,1));
      end
      fclose(fid02);

      %output slip distribution in XY for GMT
      text1=append('results/BF_',name_mod,'.sd');
      fid02=fopen(text1,'w');
      for ind=1:dimT
          %Estimate seismic moment
          smom(ind,1)=rigid*totalslip(ind,1)*area1(1)*1000000.0;
          smomd(ind,1)=rigid*totalslip(ind,1)*sind(rake(ind,1))*area1(1)*1000000.0;
          smoms(ind,1)=rigid*totalslip(ind,1)*cosd(rake(ind,1))*area1(1)*1000000.0;
          %Output to file
          fprintf(fid02,'> -Z%f\n',totalslip(ind,1));
          fprintf(fid02,'%8.4f %8.4f %8.4f %8.4f %8.4e %8.4f \n',...
              Rot2(trir(ind,1),1),Rot3(trir(ind,1),1),0.000,...
              totalslip(ind,1),smom(ind,1),rake(ind,1));
          fprintf(fid02,'%8.4f %8.4f %8.4f %8.4f %8.4e %8.4f \n',...
              Rot2(trir(ind,2),1),Rot3(trir(ind,2),1),0.000,...
              totalslip(ind,1),smom(ind,1),rake(ind,1));
          fprintf(fid02,'%8.4f %8.4f %8.4f %8.4f %8.4e %8.4f \n',...
              Rot2(trir(ind,3),1),Rot3(trir(ind,3),1),0.000,...
              totalslip(ind,1),smom(ind,1),rake(ind,1));
          fprintf(fid02,'%8.4f %8.4f %8.4f %8.4f %8.4e %8.4f \n',...
              Rot2(trir(ind,1),1),Rot3(trir(ind,1),1),0.000,...
              totalslip(ind,1),smom(ind,1),rake(ind,1));
      end
      fclose(fid02);

      %seismic moment
      smondT=sum(smomd,"all");
      smonsT=sum(smoms,"all");
      smomV=sqrt((smondT^2)+(smonsT^2));

      swV=(log10(smomV)-9.1)/1.5;
      fprintf('Mw = %8.5f \n',swV);
      fprintf('Seismic moment (Nm) = %11.5e \n',smomV);

      %Slip vectors located at centroids and seismic moment
      slip1=totalslip;
      text1=append('results/slip_BF_vectors_',name_mod,'.sd');
      fid02=fopen(text1,'w');
      for ind=1:dimT
          strike_slip1(ind,1)=totalslip(ind,1)*cosd(rake(ind,1));
          dip_slip1(ind,1)=totalslip(ind,1)*sind(rake(ind,1));
          fprintf(fid02,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4e \n',...
              ICrx(ind,1),ICry(ind,1),strike_slip1(ind,1),...
              dip_slip1(ind,1),slip1(ind,1),rake(ind,1),smom(ind,1));
      end
      fclose(fid02);

      %output displacement: observations vs calculations

      text1=append('results/Displacement_',name_mod,'.sd');
      fid02=fopen(text1,'w');
      fprintf(fid02,'Alpha= %8.4f , sa1=  %8.4f \n',alpha,sa1);
      fprintf(fid02,'Lon       Lat       EW(m)    EWe(m)   NS(m)    NSe(m)   UD(m)    UDe(m)   Mew(m)   Mns(m)   Mud(m)   rew(m)   rns(m)   rud(m)\n');
      for ind=1:1:length(obs_vy_r)
          if(horver==0 && ind<=nz0)
              fprintf(fid02,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n',...
                  stlon(ind,1),stlat(ind,1),obs_vy_r(ind,1),...
                  err_vy_r(ind,1),obs_vx_r(ind,1),err_vx_r(ind,1),...
                  -1.0*obs_vz_r(ind,1),err_vz_r(ind,1),cal_vy(ind,1),...
                  cal_vx(ind,1),0.0,Resi_vy(ind,1),Resi_vx(ind,1),0.0);
          elseif(horver==0 && ind>nz0)
              fprintf(fid02,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n',...
                  stlon(ind,1),stlat(ind,1),obs_vy_r(ind,1),...
                  err_vy_r(ind,1),obs_vx_r(ind,1),err_vx_r(ind,1),...
                  0.0,0.0,cal_vy(ind,1),cal_vx(ind,1),...
                  0.0,Resi_vy(ind,1),Resi_vx(ind,1),0.0);
          elseif(horver==1 && ind<=nz0)
              fprintf(fid02,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n',...
                  stlon(ind,1),stlat(ind,1),obs_vy_r(ind,1),...
                  err_vy_r(ind,1),obs_vx_r(ind,1),err_vx_r(ind,1),...
                  -1.0*obs_vz_r(ind,1),err_vz_r(ind,1),cal_vy(ind,1),...
                  cal_vx(ind,1),-1.0*cal_vz(ind,1),Resi_vy(ind,1),...
                  Resi_vx(ind,1),-1.0*Resi_vz(ind,1));
          elseif(horver==1 && ind>nz0)
              fprintf(fid02,'%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n',...
                  stlon(ind,1),stlat(ind,1),obs_vy_r(ind,1),...
                  err_vy_r(ind,1),obs_vx_r(ind,1),err_vx_r(ind,1),...
                  0.0,0.0,cal_vy(ind,1),...
                  cal_vx(ind,1),0.0,Resi_vy(ind,1),...
                  Resi_vx(ind,1),0.0);

          end
      end
      fclose(fid02);
    

      save strs_inv_b01
