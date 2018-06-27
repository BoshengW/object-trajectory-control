%% generate a circle trajectory then segment into n part 
function [traj,orien] = traj_circle(x,y,segment,radius)
traj = [];
points = segment;%% number of point you want to seperate your trajectory
r = radius;
beta = [];%% angle 
orien = [];%% angle 

%% get the orien,start & end point of each traj
for i = 1:points+1
    beta(i) = 360/segment*(i-1);
    target_xy = [r*cos(beta(i)/360*2*pi),r*sin(beta(i)/360*2*pi),];
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