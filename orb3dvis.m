% ***** Orb3Dvis *****

% This program is a Matlab script file written to visualise the
% output from Orbit3D which consists of the trajectory of the
% VentrAssist impeller, the pressure distribution over the surfaces
% of the conical journal bearing and the thrust bearing and a
% frequency analysis of the trajectory. The input files for this 
% script file are "plottraj.m", "plotcjb.m", and "plottb.m", which
% all must be located in the same directory as Orb3Dvis.

flaga=1;
while flaga ~= 0,
   choice1 = menu('Select option','Bearing trajectory - Eccentricity',...
      'Bearing trajectory - Rotation','CJB pressure distribution',...
      'TB pressure distribution','FFT analysis of trajectory','Quit');
   
   switch choice1
      
   case 1
      flag1 = figflag('Bearing Trajectory - Eccentricity');
      if flag1 == 0,
         h1=figure('name','Bearing Trajectory - Eccentricity');
      end
      clf(h1); plottraj;
      plot3(erx(1),ery(1),erz(1),'b+');
      hold on;
      xlabel('\epsilon_X''');
      ylabel('\epsilon_Y''');
      zlabel('\epsilon_Z');
      axis([-1.25 1.25 -1.25 1.25 -1 1]);
      grid on;
      box on;
      set(gca,'PlotBoxAspectRatio',[1.25 1.25 1])
      view(25,30);
      comet3(erx,ery,erz);
      plot3(erx,ery,erz,'LineWidth',1);
      rotate3d on;
      
   case 2
      flag2 = figflag('Bearing Trajectory - Rotation');
      if flag2 == 0,
         h2=figure('name','Bearing Trajectory - Rotation');
      end
      clf(h2); plottraj;
      plot(xrot(1),yrot(1),'b+');
      hold on;
      xlabel('\gamma_X^*');
      ylabel('\gamma_Y^*');
      axis([-1 1 -1 1]);
      box on;
      comet(xrot,yrot);
      plot(xrot,yrot,'LineWidth',1);
      
   case 3
      flag3 = figflag('CJB Pressure Distribution');
      if flag3 == 0,
         h3=figure('name','CJB Pressure Distribution');
      end     
      plotcjb;
      
   case 4
      flag4 = figflag('TB Pressure Distribution');
      if flag4 == 0,
         h4=figure('name','TB Pressure Distribution');
      end   
      plottb;
      
   case 5
      flagb = 1;
      while flagb ~= 0,
         choice2 = menu('Select option','Eccentricity ratio in X-direction',...
            'Eccentricity ratio in Y-direction','Eccentricity ratio in Z-direction',...
            'Rotation about the X-axis','Rotation about the Y-axis','Quit');
      
         switch choice2
         case 1
            flag5a = figflag('Non-dim. freq. spectrum of eccentricity in X-direction');
            if flag5a == 0,
               h5a=figure('name','Non-dim. freq. spectrum of eccentricity in X-direction');
            end
            plottraj;
            [P,F]=spectrum(erx,[],[],[],sf);
            plot(F,P(:,1));
            xlabel('Non-dimensional frequency');
            ylabel('PSD');
            axis([0 16 0 Inf]);
       
         case 2
            flag5b = figflag('Non-dim. freq. spectrum of eccentricity in Y-direction');
            if flag5b == 0,
               h5b=figure('name','Non-dim. freq. spectrum of eccentricity in Y-direction');
            end
            plottraj;
            [P,F]=spectrum(ery,[],[],[],sf);
            plot(F,P(:,1));
            xlabel('Non-dimensional frequency');
            ylabel('PSD');
            axis([0 16 0 Inf]);
         
         case 3
            flag5c = figflag('Non-dim. freq. spectrum of eccentricity in Z-direction');
            if flag5c == 0,
               h5c=figure('name','Non-dim. freq. spectrum of eccentricity in Z-direction');
            end
            plottraj;
            [P,F]=spectrum(erz,[],[],[],sf);
            plot(F,P(:,1));
            xlabel('Non-dimensional frequency');
            ylabel('PSD');
            axis([0 16 0 Inf]);
            
         case 4
            flag5d = figflag('Non-dim. freq. spectrum of rotation about X-axis');
            if flag5d == 0,
               h5d=figure('name','Non-dim. freq. spectrum of rotation about X-axis');
            end
            plottraj;
            [P,F]=spectrum(xrot,[],[],[],sf);
            plot(F,P(:,1));
            xlabel('Non-dimensional frequency');
            ylabel('PSD');
            axis([0 16 0 Inf]);
            
         case 5
            flag5e = figflag('Non-dim. freq. spectrum of rotation about Y-axis');
            if flag5e == 0,
               h5e=figure('name','Non-dim. freq. spectrum of rotation about Y-axis');
            end
            plottraj;
            [P,F]=spectrum(yrot,[],[],[],sf);
            plot(F,P(:,1));
            xlabel('Non-dimensional frequency');
            ylabel('PSD');
            axis([0 16 0 Inf]);   
            
         case 6
            flag5a = figflag('Non-dim. freq. spectrum of eccentricity in X-direction');
            if flag5a == 1,
               close(h5a);
            end
            flag5b = figflag('Non-dim. freq. spectrum of eccentricity in Y-direction');
            if flag5b == 1,
               close(h5b);
            end
            flag5c = figflag('Non-dim. freq. spectrum of eccentricity in Z-direction');
            if flag5c == 1,
               close(h5c);
            end
            flag5d = figflag('Non-dim. freq. spectrum of rotation about X-axis');
            if flag5d == 1,
               close(h5d);
            end
            flag5e = figflag('Non-dim. freq. spectrum of rotation about Y-axis');
            if flag5e == 1,
               close(h5e);
            end            
            flagb = 0;
         
         end     
         
      end
     
   case 6    
      flag1 = figflag('Bearing Trajectory - Eccentricity');
      if flag1 == 1,
         close(h1);
      end
      flag2 = figflag('Bearing Trajectory - Rotation');
      if flag2 == 1,
         close(h2);
      end      
      flag3 = figflag('CJB Pressure Distribution');
      if flag3 == 1,
         close(h3);
      end
      flag4 = figflag('TB Pressure Distribution');
      if flag4 == 1,
         close(h4);
      end        
      flaga = 0;         
      
   end
   
end

