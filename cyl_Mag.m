function [Bz,Bp,Gradient,M_object,F_object] = cyl_Mag(x0,y0,D,L,orien_mag,i)
    syms z p;
    %%constant in magentic field model
    C = 1.45/(4*pi); %% Br/4pi based on N52: Br = 1.45T
    a = L/2;  %% radius of permanent magnet
    slope = tan(orien_mag/360*2*pi);  %% convert orientation into linear function 
    %% position of external magnet 
    x1 = D*cos(orien_mag/360*2*pi);
    y1 = D*sin(orien_mag/360*2*pi); 
    %% magnetic field model 
    B_z = @(R,phi,z,p)(R*(L/2-z)./(R.^2+(L/2-z)^2+p^2-2*R*p.*cos(phi)).^1.5 + (R*(L/2+z))./(R.^2+(L/2+z)^2+p^2-2*R*p.*cos(phi)).^1.5);
    B_p = @(R,phi,z,p)(-R.*(2*p-2*R.*cos(phi))/(2*(R.^2+(L/2-z)^2+p^2-2*R*p.*cos(phi)).^1.5)+R.*(2*p-2*R.*cos(phi))./(2*(R.^2+(L/2+z)^2+p^2-2*R*p.*cos(phi)).^1.5));
    %% magnetic gradient 
    dz_Bz = matlabFunction(diff(B_z,z,1));
    dp_Bz = matlabFunction(diff(B_z,p,1));
    dz_Bp = matlabFunction(diff(B_p,z,1));
    dp_Bp = matlabFunction(diff(B_p,p,1));
    %% axial and radial distance 
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
            Bz = @(R,phi)(B_z(R,phi,z_m1,p_m1));
            Bp = @(R,phi)(B_p(R,phi,z_m1,p_m1));
            field_z = -C*dblquad(Bz,0,a,0,2*pi);
            field_p = -C*dblquad(Bp,0,a,0,2*pi);
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
            Bz = @(R,phi)(B_z(R,phi,z_m1,p_m1));
            Bp = @(R,phi)(B_p(R,phi,z_m1,p_m1));
            field_z = -C*dblquad(Bz,0,a,0,2*pi);
            field_p = -C*dblquad(Bp,0,a,0,2*pi);
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

    %% calculation of gradient of the below magnet 
        Bz_dz = @(R,phi)(dz_Bz(R,p_m1,phi,z_m1));
        Bz_dp = @(R,phi)(dp_Bz(R,p_m1,phi,z_m1)); 
        Bp_dz = @(R,phi)(dz_Bp(R,p_m1,phi,z_m1));
        Bp_dp = @(R,phi)(dp_Bp(R,p_m1,phi,z_m1));

        gradient_z_m1_z = -C*dblquad(Bz_dz,0,a,0,2*pi);
        gradient_z_m1_p = -C*dblquad(Bz_dp,0,a,0,2*pi);
        gradient_p_m1_z = -C*dblquad(Bp_dz,0,a,0,2*pi);
        gradient_p_m1_p = -C*dblquad(Bp_dp,0,a,0,2*pi);   
        Gradient = [gradient_z_m1_z,gradient_p_m1_z;gradient_z_m1_p,gradient_p_m1_p];
       %% Magnetic force of object 
        %% dipole moment: M_object = M*V 
        %% Magnitization M = Br/u0: Br = 1.32T in N42; permeability of vacuum u0 = 4*pi*10^-7 H/m
        %% Volumn of object V = 4/3*pi*R^3
        M_object =  1.32/(4*pi)*10^7*4/3*pi*(1.5*10^-3)^3*M_object;
         %% real magnitization of object
        F_object_zp = (M_object*Gradient); %% magnetic force in z&p space [Fz,Fp]
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
     
        F_object = [F_object_x,F_object_y];
end