function [traj,orien] = traj_S(x,y,segment,radius)
%%% this fucntion used for generate a "S" traj
traj = [];
points = segment;%% number of point you want to seperate the traj
r = radius;
beta = [];%% angle 
orien = [];%% angle 
%% centroid 1
cx1 = -0.6*10^-2; 
cy1 = 0;
%% centroid 2
cx2 = 0.6*10^-2;
cy2 = 0;
%% generate the traj
for i = 1:points+1
    beta(i) = 360/segment*(i-1);
    if i<=(points+1)/2
        target_xy = [cx1+r*cosd(180-beta(i)),cy1+r*sind(180-beta(i))];
    else
        target_xy = [cx2+r*cosd(beta(i)),cy2+r*sind(beta(i))];
    end
    
    traj = [traj;target_xy];
    if i>=2       
        traj_orien = (traj(i,:) - traj(i-1,:))/norm(traj(i,:) - traj(i-1,:));
        traj_orien = atan2d(traj_orien(2),traj_orien(1));
        if traj_orien<0&&traj_orien>=-90 
            traj_orien = 360 - abs(traj_orien);       
        end
        if traj_orien<-90&&traj_orien>=-180 
            traj_orien = 360- abs(traj_orien);
        end
        orien = [orien;traj_orien];        
    end  
end
traj_orien = (traj(1,:) - [x,y])/norm(traj(1,:) - [x,y]);
traj_orien = atan2d(traj_orien(2),traj_orien(1));
orien = [traj_orien;orien];
target_xy = [x,y];
traj = [target_xy;traj];


end