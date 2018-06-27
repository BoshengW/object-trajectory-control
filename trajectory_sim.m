%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   trajectory control simulation
%   Bosheng Wu bw493@drexel.edu
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms z p;
%% magnetic field model
%% parameter set
L = 2.54*10^-2;   %% length of below magnet
a = L/2;          %% radius of below magnet 
Br = 1.45/(4*pi); %% magnitization of spherical magnet
w_magnet = 0;     %% intial rotation speed of magnet    
orien_mag = 0;        %% initial orientation of magnet
slope = tan(orien_mag/360*2*pi);  %% convert orientation into linear function 
K = 1000;         %% number of iteration
T = 1;          %% time of each iteration (sec)
D = 4.6*10^-2;      %% raduis from magnet to center of coordinate
r = 3*10^-3/2;    %% radius of the object
error2 = [];

%% test 
V = [];
distance =  [];
F = [];
%%
%% parameter of control loop
error = zeros(K);
Kp = 3;   %% P
Kd = 3;       %% D
info_angle = [];
info_angle = [info_angle,orien_mag]; %% initialize

%% orientation of object
theta_ob = 0;
w_object = 0;
%% initial position of object 
    %% object position(centroid)
    x0 = -1.2*10^-2;
    y0 = -0.1*10^-2;

    %% below magnet position(centroid)
    x1 = D*cos(orien_mag/360*2*pi);
    y1 = D*sin(orien_mag/360*2*pi); 

%% magnetic field model
%% magnetic field 
B_z = @(R,phi,z,p)(R*(L/2-z)./(R.^2+(L/2-z)^2+p^2-2*R*p.*cos(phi)).^1.5 + (R*(L/2+z))./(R.^2+(L/2+z)^2+p^2-2*R*p.*cos(phi)).^1.5);
B_p = @(R,phi,z,p)(-R.*(2*p-2*R.*cos(phi))/(2*(R.^2+(L/2-z)^2+p^2-2*R*p.*cos(phi)).^1.5)+R.*(2*p-2*R.*cos(phi))./(2*(R.^2+(L/2+z)^2+p^2-2*R*p.*cos(phi)).^1.5));
%% magnetic gradient 
dz_Bz = matlabFunction(diff(B_z,z,1));
dp_Bz = matlabFunction(diff(B_z,p,1));
dz_Bp = matlabFunction(diff(B_p,z,1));
dp_Bp = matlabFunction(diff(B_p,p,1));

%% generate desired trajectory 
[traj,orien] = traj_S(x0,y0,36,0.6*10^-2);
% [traj,orien] = traj_circle(x0,y0,36,2*10^-2);
segment = size(orien,1);
step = 1;
check = true;
limit = 0.1*10^-2;
cross_error =0;
%% plot desired trajectory
group = [];
group1 = [];
figure(1)
plot(traj(:,1),traj(:,2),'b.')

xlabel("X axis: (m)")
ylabel("Y axis: (m)")
hold on
%% iteration begin
i =1;
distance = [];
%% Cycle for the simulation
while(step<38)%(step<37)
   
   %% code plot the trajectory
   plot(traj(:,1),traj(:,2),'b.')
   hold on
   plot(x0,y0,'r.')
   axis([-0.03,0.03,-0.03,0.03])
   hold on
%    plot(x1,y1,'bo')
%    hold on

%%********************* fucntion cyl_Mag***************%%
   %% distance from object to side magnet
    %% check the location of object 
    p_m1 = abs(slope*x0-y0)/sqrt(slope^2+1);
    z_m1 = sqrt((D*cos(orien_mag/360*2*pi)-x0)^2+(D*sin(orien_mag/360*2*pi)-y0)^2-p_m1^2);
    %% check the sign of p_m1
    if x1<0
        if slope*x0-y0>=0
            p_m1 = p_m1;
        else
            p_m1 = -p_m1;
        end
    else if x1>0
          if slope*x0-y0>=0
            p_m1 = -p_m1;
          else
            p_m1 = p_m1;
          end 
        end
    end
    if x1==0 && y1>0
        if x0<=0
            p_m1 = p_m1;
        else 
            p_m1 = -p_m1;
        end        
    end
    if x1==0 && y1<0
        if x0<=0
            p_m1 = -p_m1;
        else 
            p_m1 = p_m1;
        end        
    end
      
    %% magnetization of the object 
        %% initial position of obeject 
        if i==1
            object_z = @(R,phi)(B_z(R,phi,z_m1,p_m1));
            object_p = @(R,phi)(B_p(R,phi,z_m1,p_m1));
            field_z = -Br*dblquad(object_z,0,a,0,2*pi);
            field_p = -Br*dblquad(object_p,0,a,0,2*pi);
            B_object = [field_z,field_p];
                   
        %% initial the magnitization into xy
        %% convert the magnetization from zp space into xy space
            if p_m1>0 && x1>0
                B_x = -field_z*cos(orien_mag/360*2*pi) - field_p*sin(orien_mag/360*2*pi);
                B_y = -field_z*sin(orien_mag/360*2*pi) + field_p*cos(orien_mag/360*2*pi);
            end
            %%if it is negative, x1>0
            if p_m1<0 && x1>0
                B_x = -field_z*cos(orien_mag/360*2*pi) - field_p*sin(orien_mag/360*2*pi);
                B_y = -field_z*sin(orien_mag/360*2*pi) + field_p*cos(orien_mag/360*2*pi);
            end
            %%if it is positive, x1<0
            if p_m1>0 && x1<0
                B_x = -field_z*cos(orien_mag/360*2*pi) + field_p*sin(orien_mag/360*2*pi);
                B_y = -field_z*sin(orien_mag/360*2*pi) - field_p*cos(orien_mag/360*2*pi);
            end
            %%if it is negative, x1<0
            if p_m1<0 && x1<0
                B_x = -field_z*cos(orien_mag/360*2*pi) - field_p*sin(orien_mag/360*2*pi);
                B_y = -field_z*sin(orien_mag/360*2*pi) + field_p*cos(orien_mag/360*2*pi);
            end
            if p_m1 ==0
                B_x = -field_z*cos(orien_mag/360*2*pi) + field_p*sin(orien_mag/360*2*pi);
                B_y = -field_z*sin(orien_mag/360*2*pi) + field_p*cos(orien_mag/360*2*pi);
            end
            B_object_xy = [B_x,B_y];
            B_xy = B_object_xy/norm(B_object_xy);
            B_angle = atan2d(B_xy(2),B_xy(1));
            M_object_xy = B_xy;
            M_angle = atan2d(B_y,B_x);
            if M_angle<0&&M_angle>=-90 
            M_angle = 360 - abs(M_angle); 
            end
            if M_angle<-90&&M_angle>=-180 
            M_angle = 360- abs(M_angle);
            end 
        %% After initialize position of object 
        else 
            object_z = @(R,phi)(B_z(R,phi,z_m1,p_m1));
            object_p = @(R,phi)(B_p(R,phi,z_m1,p_m1));
            field_z = -Br*dblquad(object_z,0,a,0,2*pi);
            field_p = -Br*dblquad(object_p,0,a,0,2*pi);
            B_object = [field_z,field_p];
            
            %% convert Bz,Bp into xy
            if p_m1>0 && x1>0
                B_x = -field_z*cos(orien_mag/360*2*pi) - field_p*sin(orien_mag/360*2*pi);
                B_y = -field_z*sin(orien_mag/360*2*pi) + field_p*cos(orien_mag/360*2*pi);
            end
            %%if it is negative, x1>0
            if p_m1<0 && x1>0
                B_x = -field_z*cos(orien_mag/360*2*pi) - field_p*sin(orien_mag/360*2*pi);
                B_y = -field_z*sin(orien_mag/360*2*pi) + field_p*cos(orien_mag/360*2*pi);
            end
            %%if it is positive, x1<0
            if p_m1>0 && x1<0
                B_x = -field_z*cos(orien_mag/360*2*pi) + field_p*sin(orien_mag/360*2*pi);
                B_y = -field_z*sin(orien_mag/360*2*pi) - field_p*cos(orien_mag/360*2*pi);
            end
            %%if it is negative, x1<0
            if p_m1<0 && x1<0
                B_x = -field_z*cos(orien_mag/360*2*pi) - field_p*sin(orien_mag/360*2*pi);
                B_y = -field_z*sin(orien_mag/360*2*pi) + field_p*cos(orien_mag/360*2*pi);
            end
            if p_m1 ==0
                B_x = -field_z*cos(orien_mag/360*2*pi) + field_p*sin(orien_mag/360*2*pi);
                B_y = -field_z*sin(orien_mag/360*2*pi) + field_p*cos(orien_mag/360*2*pi);
            end
            B_object_xy = [B_x,B_y];
            B_xy = B_object_xy/norm(B_object_xy);
            B_angle = atan2d(B_xy(2),B_xy(1));
        end
            %% adjust the value of B_angle 
            if B_angle<0&&B_angle>=-90 
            B_angle = 360 - abs(B_angle); 
            end
            if B_angle<-90&&B_angle>=-180 
            B_angle = 360- abs(B_angle);
            end
        %% magnetization of the object
        if i==1
            M_object = B_object/norm(B_object);
        else
            intsect = M_angle-B_angle; %% angle between new B and old M of object
            M_angle = M_angle + rotate_angle; %% the object angle in xy space 
            if M_angle>360
                M_angle = M_angle-360;
            else
                if M_angle<0
                    M_angle = 360 + M_angle;
                end
            end
            
            if abs(rotate_angle)>=abs(intsect)
                M_angle = B_angle;
                Bzp_angle = atan2d(B_object(1),B_object(2));      
            else
                Bzp_angle = atan2d(B_object(1),B_object(2));               
            end
            %% adjust the angle of Bzp_angle
            if Bzp_angle<0&&Bzp_angle>=-90 
            Bzp_angle = 360 - abs(Bzp_angle); 
            end
            if Bzp_angle<-90&&Bzp_angle>=-180 
            Bzp_angle = 360- abs(Bzp_angle);
            end
            M_object = [sin((Bzp_angle+M_angle-B_angle)/360*2*pi),cos((Bzp_angle+M_angle-B_angle)/360*2*pi)];
        end
%         %%%%%%%%%%%%%%%%%%%%%test if no torque
         M_object = B_object/norm(B_object);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%
%          display([B_angle,M_angle]);
    %% calculation of gradient of the below magnet 
        Bz_dz = @(R,phi)(dz_Bz(R,p_m1,phi,z_m1));
        Bz_dp = @(R,phi)(dp_Bz(R,p_m1,phi,z_m1)); 
        Bp_dz = @(R,phi)(dz_Bp(R,p_m1,phi,z_m1));
        Bp_dp = @(R,phi)(dp_Bp(R,p_m1,phi,z_m1));

        gradient_z_m1_z = -Br*dblquad(Bz_dz,0,a,0,2*pi);
        gradient_z_m1_p = -Br*dblquad(Bz_dp,0,a,0,2*pi);
        gradient_p_m1_z = -Br*dblquad(Bp_dz,0,a,0,2*pi);
        gradient_p_m1_p = -Br*dblquad(Bp_dp,0,a,0,2*pi);   
        gradient_object = [gradient_z_m1_z,gradient_p_m1_z;gradient_z_m1_p,gradient_p_m1_p];
    %% Magnetic force of object 
    display(M_object);
    M_object = 0.015*M_object;
     %% real magnitization of object
    F_object_zp = (M_object*gradient_object); %% magnetic force in z&p space [Fz,Fp]
    %%if p_m1 is positive, x1>0 
    if p_m1>0 && x1>0
        F_object_x = -F_object_zp(1)*cos(orien_mag/360*2*pi) + F_object_zp(2)*sin(orien_mag/360*2*pi);
        F_object_y = -F_object_zp(1)*sin(orien_mag/360*2*pi) + F_object_zp(2)*cos(orien_mag/360*2*pi);
    end
    %%if it is negative, x1>0
    if p_m1<0 && x1>0
        F_object_x = -F_object_zp(1)*cos(orien_mag/360*2*pi) - F_object_zp(2)*sin(orien_mag/360*2*pi);
        F_object_y = -F_object_zp(1)*sin(orien_mag/360*2*pi) + F_object_zp(2)*cos(orien_mag/360*2*pi);
    end
    %%if it is positive, x1<0
    if p_m1>0 && x1<0
        F_object_x = -F_object_zp(1)*cos(orien_mag/360*2*pi) + F_object_zp(2)*sin(orien_mag/360*2*pi);
        F_object_y = -F_object_zp(1)*sin(orien_mag/360*2*pi) - F_object_zp(2)*cos(orien_mag/360*2*pi);
    end
    %%if it is negative, x1<0
    if p_m1<0 && x1<0
        F_object_x = -F_object_zp(1)*cos(orien_mag/360*2*pi) - F_object_zp(2)*sin(orien_mag/360*2*pi);
        F_object_y = -F_object_zp(1)*sin(orien_mag/360*2*pi) + F_object_zp(2)*cos(orien_mag/360*2*pi);
    end
    if p_m1 ==0
        F_object_x = -F_object_zp(1)*cos(orien_mag/360*2*pi) + F_object_zp(2)*sin(orien_mag/360*2*pi);
        F_object_y = -F_object_zp(1)*sin(orien_mag/360*2*pi) + F_object_zp(2)*cos(orien_mag/360*2*pi);
    end
    
%%*********************end of  fucntion cyl_Mag***************%%


    %% Magnetic torque of the object 
    B_3D = [B_object_xy,0];
    M_3D = [M_object_xy,0];
    torque = cross(M_3D,B_3D); %% get the magnitude of torque on object 
     %% angular velocity 
     w = sign(torque(3))*norm(torque)/0.002;
    %% update object state 
    %% (x',y')speed of object in x&y (refer to model based on linear motion experiment data)
% 
        v_x = sign(F_object_x)*(abs(F_object_x)-0.005)/0.14*10^-3;  
        v_y = sign(F_object_y)*(abs(F_object_y)-0.005)/0.14*10^-3;

        if abs(F_object_x)<=0.005
            v_x = 0;
        end
        if abs(F_object_y)<=0.005
            v_y = 0;
        end
    %% update object ouput
     %% position of object 
        x0 = x0 + v_x*T;
        y0 = y0 + v_y*T;
        rotate_angle = w*57.3*T; %% convert rad to degree 1 rad = 57.3 degree
        [orien_mag,step,check,cross_error] = traj_carrot(x0,y0,step,traj,orien,limit,check,segment,Kp,Kd,cross_error);  
        display(orien_mag);
        pause(0.5);
    %% update the magnet position
    x1 = D*cos(orien_mag/360*2*pi);
    y1 = D*sin(orien_mag/360*2*pi); 
    if step>2
        limit = 0.05*10^-2;
    end
    slope = tan(orien_mag/360*2*pi);
  
    i= i+1;

    
end
