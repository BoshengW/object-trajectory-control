function [orien_mag,step,check,cross_error] = traj_carrot(x,y,step,traj,orien,limit,check,segment,Kp,Kd,cross_error)
  
  part = step;
  distance = sqrt((x-traj(part+1,1))^2+(y-traj(part+1,2))^2);
  %% check if the object reach the final points in this 
    %% if reach then switch to next part  
      if distance<= limit
          part = part+1;%% switch to next part
          step = part;
          check = true;
      end
      if part>segment
          return;
      end
  slope_L = tand(orien(part)); %% slope of this part 
  b = traj(part,2) - slope_L*traj(part,1);
  
  %% get the error between real & desired trajectory
  if check
      check = false;
      error_rate = 0;
      

  else
      last_error = cross_error;
      if step<20
          cross_error = sign(slope_L*x-y+b)*abs(slope_L*x-y+b)/sqrt(1+slope_L^2);
          %cross_error = -sign(slope_L*x-y+b)*abs(slope_L*x-y+b)/sqrt(1+slope_L^2); %% track cross error
      else
         cross_error = sign(slope_L*x-y+b)*abs(slope_L*x-y+b)/sqrt(1+slope_L^2); %% track cross error 
      end
    
 

  display(atand(cross_error/distance));
  end
  %% proportional control 
 orien_mag = orien(part) + Kp*atand(cross_error/distance);
 
%     orien_mag =  orien(part) + (Kp*cross_error)*57.3;
    
    %% following part is to limit the max orientation compenstation
%     max = orien(part) + 90;
%     min = orien(part) - 90;
%     if orien_mag<min
%         orien_mag = min;
%     end
%     if orien_mag>max
%         orien_mag = max;
%     end
%     if orien_mag <0
%       orien_mag = 360 + orien_mag;
%     end

end