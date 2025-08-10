      clear; close all;
      %A. Meneses-Gutierrez & T. Saito (2025)
      %--------------------------------------------------------------------------
      % This program computes the slip distribution corresponding to the 
      % basis function (basisf) in the given fault geometry.

      % This code was developed to produce results presented in:
      % "Linking Coseismic Slip and Afterslip in Intraplate Earthquakes:
      % A Case Study of the 2016 Central Tottori Earthquake, Japan"
      % (2025JB031677).

      %Stress drop basis functions
 

      load('faultP_par01','disy','disx','disz','trir1','trir2','trir3',...
          'dimT','dimk');
      load('area_pos_B01','basisf');%Basis functions distribution
      load('surface_disp_T01','ux_slip','uy_slip','uz_slip');
      load('meshresponsephs_T01','slip_traction');
      
      tic

      size1=size(basisf);
      num_base=size1(2);

      G=slip_traction; %MPa/m
      invG=inv(G); %m/MPa

      for ind = 1:num_base
        d1=basisf(:,ind); %MPa
        mm=invG*d1;
        slipG_f0(:,ind)=mm; % m 

        ux_basis(:,ind)=ux_slip*slipG_f0(:,ind);
        uy_basis(:,ind)=uy_slip*slipG_f0(:,ind);
        uz_basis(:,ind)=uz_slip*slipG_f0(:,ind);

      end
      
      caltime=toc

      save stressbase_T01
