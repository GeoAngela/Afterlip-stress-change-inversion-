      clear; close all;
      %A. Meneses-Gutierrez & T. Saito (2025)
      %--------------------------------------------------------------------------
      % This program estimates the traction on each subfault due to 1 m of slip
      % in the specified slip direction.

      % This program uses the function **TDstressHS** from the software
      % developed by Nikkhoo and Walter (2015) to compute strains and
      % stresses due to triangular dislocations in a half-space.

      % This code was developed to produce results presented in:
      % "Linking Coseismic Slip and Afterslip in Intraplate Earthquakes:
      % A Case Study of the 2016 Central Tottori Earthquake, Japan"
      % (2025JB031677).

      %Z is positive downward

      tic

      load('slip_vector_T01','vx','vy','vz');
      load('faultP_par01','disy','disx','disz','trir1','trir2','trir3',...
          'dimT','dimk');
      load('faultP_par01','X','Y','Z');
      load('slip_vector_T01','n');
      load('slip_vector_T01','rake');%slip direction for calculation

      %Elastic constants
      rigid=30.0;
      poi=0.25;
      Kel=2*rigid*(1+poi)/(3*(1-2*poi));
      miu=3*Kel*(1-2*poi)/(2*(1+poi));
      lambda=3*Kel*poi/(1+poi);

  %%  Calculate stress distribution from each mesh
      trir=[trir1(1:dimT,1), trir2(1:dimT,1),trir3(1:dimT,1)];
      dis=[disy(1:dimk,1), disx(1:dimk,1), disz(1:dimk,1)];

      for i=1:1:size(trir,1)
          % Note the coordinates are different
          % Stress:
          % Calculated stress tensor components in EFCS. The six columns of Stress
          % are Sxx, Syy, Szz, Sxy, Sxz and Syz, respectively. The stress components
          % have the same unit as Lame constants.
          %% m/km*GPa = MPa
          [stress,strain] = ...
              TDstressHS(Y,X,-1*Z,[dis(trir(i,1),1) dis(trir(i,1),2) -1*dis(trir(i,1),3)],...
              [dis(trir(i,2),1) dis(trir(i,2),2) -1*dis(trir(i,2),3)],...
              [dis(trir(i,3),1) dis(trir(i,3),2) -1*dis(trir(i,3),3)],...
              cosd(rake(i,1)),sind(rake(i,1)),0,miu,lambda);

       if(mod(i,100) == 0)
        i
        toc
       end

      S_xx=stress(:,2);
      S_xy=stress(:,4);
      S_xz=-1*stress(:,6);
      S_yy=stress(:,1);
      S_yz=-1*stress(:,5);
      S_zz=stress(:,3);

      val2 = vx.*(n(:,1).*S_xx + n(:,2).*S_xy + n(:,3).*S_xz) + ...
               vy.*(n(:,1).*S_xy + n(:,2).*S_yy + n(:,3).*S_yz) + ...
               vz.*(n(:,1).*S_xz + n(:,2).*S_yz + n(:,3).*S_zz);
      slip_traction(:,i)=val2'; % MPa/m

      end %% for each slip on each mesh


      caltime = toc
      save meshresponsephs_T01

