function [f_s,NT_Level,p_p,a_s,r_r,NT_Q,NT_depth,NT_vel,CN,output] =BRAHMA_DRFF_func(inflow,advanced_spill,z,st,n1)

d_s=3.28;                                                                   %dead storage in MCM
g_s=7.71;                                                                   %gross storage in MCM
F_P=135;                                                                    %firm power in MW
TWL=246.45;                                                                 %Tail water level in m
eeta=0.8;                                                                   %plant_efficiency
i_s(1,1)=7.71;                                                              % initial storage  
P=405;                                                                      % Plant capacity                             
g=9.81;                                                                     % acceleration due to gravity
delt=2;                                                                     % time step
cons=.01;                                                                   % artificial viscosity
delx=500;                                                                   % del x
T=3600;
r=(T/delt);
YY=146;                                                                    % number of grids
L=YY*delx;

val= xlsread('Bathymetry.xlsx');                                           % the bathymetry sheet is available with three columns: first column=chaingae, second column is the width and third column is the bathymetry.
[AA,EE]=size(val);

m=1;                                                                       % variable for interpolating computational section
x(1)=0;                                                                    % chainage of starting point
K(1,1)=0.03;                                                               % coefficient of contraction-expansion
n(1,1)=0.035;                                                              % manning's coefficient
Ba(1,1)= val(1,3);                                                         % bathymetry
Ba(AA,1)=val(AA,3);
R =(val(AA,1)-val(1,1))/delx;                                              % number of sections
b(1,1)= val(1,2);                                                          % first element of width of the river              
b(AA,1)=val(AA,2);                                                         % width of the river
Cc1(1,1)=val(1,1);                                                         % first element of chainage
Cc1(AA,1)=val(AA,1);                                                       % chainage

cc1_=28;
cc1_1=69;
cc1_2=146;
inflow1=inflow*3600*10^-6;                                                  % inflow to reservoir
for t=1:z
    nu(t,1)=3;                                                              % three units of power
    elev(t,1)=0.0049*i_s(t,1)^3-0.2129*i_s(t,1)^2+3.318*i_s(t,1)+549.54;    % elevation
    area(t,1)=(6*10^-5*i_s(t,1)^3)-0.0035*i_s(t,1)^2+0.0949*i_s(t,1)+0.1237;% area 
    mef(t,1)=6.9*3600*10^-6;                                                %environmental flow
    a_s(t,1)=i_s(t,1)+inflow1(t,1)-d_s;                                     %active storage
    p_h(t,1)=elev(t,1)-TWL;                                                 %power head
    p_d(t,1)=(n1(t,1)*3600*1000*nu(t,1)*F_P)/(9.8*10^6*p_h(t,1)*eeta); 
    if a_s(t,1)>p_d(t,1)
        p_r(t,1)=p_d(t,1);                                                 %power release
    else
        p_r(t,1)=a_s(t,1);
    end
    p_p(t,1)=(9.8*10^6*p_r(t,1)*p_h(t,1)*eeta)/(3600*1000);                 %power produce 
  
    if a_s(t,1)+d_s-p_r(t,1)-mef>g_s
        spill(t,1)=a_s(t,1)+d_s-p_r(t,1)-g_s;                               %spill
    else
        spill(t,1)=0;
    end
    if a_s(t,1)+d_s-p_r(t,1)-mef(t,1)>g_s
        f_s(t,1)=g_s;                                                      %final storage
    else
        f_s(t,1)=a_s(t,1)+d_s-p_r(t,1)-mef(t,1);
    end
    F_S(t,1)=f_s(t,1)-advanced_spill(t,1);

    i_s(t+1,1)=F_S(t,1);
    r_r(t,1)=(spill(t,1)+mef(t,1)+advanced_spill(t,1))/3600*10^6;           % rereservoir release
end

% % % % % % % % % % % % % % % without online
data_1=xlsread('sensor reading.xlsx',1);                                   %sensor reading of three sensors (hypothetical)
data_2=xlsread('sensor reading.xlsx',2);
data_3=xlsread('sensor reading.xlsx',3);

% Online Sensor data from cloud
% % % % % % % % % % % % % % % 
%channelid=;
%channelapi=;
%channeldata=thingSpeakRead(channelid, 'ReadKey', channelapi);
%    
%chdata=channeldata;

%%%%%%%%  without online data
 aa2=data_1(:,2);
 aa3=data_2(:,2);                              
 aa4=data_3(:,2);

%%%%%%%%  with online data
%aa2=data_1;
%aa3=chdata*10^-2;
%aa4=data_3(:,2);
%%%%%%%%%
t11_1=aa2(end)                                                              % extracting the last value of 1st sensor
t11_2=aa3(end)                                                              % extracting the last value of 2nd sensor        
t11_3=aa4(end)                                                              % extracting the last value of 3rd sensor 

S2=53;                                                                      % sensor location sections
S3=94;
S4=146;

% % % % % % % % % % % % % 

   %BRAHMA 1D HYDRODYNAMIC MODEL
     for i=2:AA
        
        for m=m:R
           x(m+1)= x(m)+delx;
            
              if x(m+1)<=val(i,1)
                                  
                  Ba(m+1,1)=val(i-1,3)+((val(i,3)-val(i-1,3))/(val(i,1)-val(i-1,1)))*(x(m+1)-val(i-1,1));   %bathymetry interpolation
                  b(m+1,1)=val(i-1,2)+((val(i,2)-val(i-1,2))/(val(i,1)-val(i-1,1)))*(x(m+1)-val(i-1,1));    % width interpolation
                  Cc1(m+1,1)=val(i-1,1)+((val(i,1)-val(i-1,1))/(val(i,1)-val(i-1,1)))*(x(m+1)-val(i-1,1));  % chainage interpolation
                   if Ba(m+1)<=Ba(m)
                  sx(m+1,1)=-(Ba(m+1)-Ba(m))/(x(m+1)-x(m));                                                  % slope
%                   K(m+1,1)=val(i,6)-((val(i,5)-val(i-1,6))/(val(i,1)-val(i-1,1)))*(val(i,1)-x(m+1));
%                   n(m+1,1)=val(i-1,4)+((val(i,4)-val(i-1,4))/(val(i,1)-val(i-1,1)))*(x(m+1)-val(i-1,1));
                   
              else
                     sx(m+1,1)=0.0001;
                   end
                   m=m+1;
              else
               break;          
              end
    
        end      
    
         
     end
      sx(1,1)=sx(2,1); 

p=1;
% % % % % % % % % % % % % % % % % 
          
            aa1=xlsread('RL of sensors.xlsx');                              % RL of sensor
            bb1(:,1)=Cc1(:,1);                                              % chainage interval 
            cc1=aa1(:,1);                                                   %cc1=chainage of sensor location
            we1=aa1(:,2);                                                   % RL of three sensors                                
            yy1=bb1(:,1);                                                   % chainage interval
            mpp1=interp1(cc1,we1,yy1,'linear','extrap');                    % interpolation of RL of all sections
            

            t_2_4=[t11_1,t11_2,t11_3];
            t_2_4=transpose(t_2_4);                                         % sensor values

            hpp1=interp1(cc1,t_2_4,yy1,'linear','extrap');                  % interpolation of sensor values 
            surfelv=mpp1(:,:);
            depth=mpp1(:,1)-Ba(:,1);
            depth1=depth(:,1)-hpp1(:,1);
            dep=depth1(:,:);
            elev=dep+Ba;

                                                                            % Assuming the initial flow of the river with one sensor(one sensor is in operation)
                                                                            % discharge from sensor depth
                q3=(elev(S3,1)/89.79)^(1/0.00759);                          % stage-discharge relationship
                q1=q3;
                q2=q3;
                q4=q3;
            
        for  i=1:S2
       q(i,1)=q1;
       h(i,1)=(q(i,1)*n(1,1)/(b(i,1)*sx(i,1)^(1/2)))^(3/5);                  %initial depth upto sensor 2
       a(i,1)=b(i,1)*h(i,1);
       u(i,1)=q(i,1)/a(i,1);
        end

        for  i=S2:S3
       q(i,1)=q2;
       h(i,1)=(q(i,1)*n(1,1)/(b(i,1)*sx(i,1)^(1/2)))^(3/5);
       h(S2,:)=dep(S2,:);                                                   %initial depth upto from sensor 2 to sensor 3 
       a(i,1)=b(i,1)*h(i,1);
       u(i,1)=q(i,1)/a(i,1);
        end
        
         for  i=S3:146
       q(i,1)=q3;
       h(i,1)=(q(i,1)*n(1,1)/(b(i,1)*sx(i,1)^(1/2)))^(3/5);
       h(S2,:)=dep(S2,:);
       h(S3,:)=dep(S3,:);
       h(S4,:)=dep(S4,:);
       a(i,1)=b(i,1)*h(i,1);
       u(i,1)=q(i,1)/a(i,1);
         end 

         xlswrite('initial_depth.xlsx',h,'depth');                          % initial depth 
         xlswrite('initial_depth.xlsx',u,'velocity');                       % initial velocity

% % % % % % % % % % % % % % % % % % % 

pp=r_r;                                                                    % reservoir release for routing 
[zx,zy]=size(pp);


for i=1:zx 
    xx(i,1)=pp(i);
end
xx(i+1,1)=pp(i);
 for i=1:zx
   for k=p:(i*r+1)
   q(1,k)=xx(i)+((xx(i+1)-xx(i))/r)*(k-p);                                  % interpolation accourding to time step
   end
   p=(i*r+1);                                                               
 end

for k=1:p
    for i=1
      h(i,k)=((q(1,k)*n(1,1)/(b(i,1)^(8/3)*sx(1,1)^(1/2)))^(3/5)*(1+.855*(q(1,k)*n(1,1)/(b(i,1)^(8/3)*sx(1,1)^(1/2)))^(3/5)))*b(i,1);
      u(i,k)= q(1,k)/(b(i,1)*h(i,k));
    end
    a(1,k)=h(1,k)*b(1,1);
    q(1,k)=a(1,k)*u(1,k);
end
                             
       new_h =xlsread('initial_depth.xlsx',2); 
       new_u =xlsread('initial_depth.xlsx',3);
  for  i=1:YY  
          h(i,1)=new_h(i,1);
        a(i,1)=b(i,1)*h(i,1);
        u(i,1)=new_u(i,1);
        q(i,1)=a(i,1)/u(i,1);
       
end
      
for i=1:YY
            hmax(i,1)=0;
            vmax(i,1)=0;
            v_hmax(i,1)=0;
end
  
 for k=1:p-1

       for i=2:24-1
            ap(i,k)=a(i,k)-(delt/(delx))*(q(i,k)-q(i-1,k));
            qp(i,k)=q(i,k)+delt*g*a(i,k)*(sx(i,1)-(n(1,1)^2*q(i,k)^2*(b(i,1)+2*h(i,k))^(4/3)*a(i,k)^(-10/3)))-(delt/delx)*(q(i,k)^2/a(i,k)+g*a(i,k)*h(i,k)/2-q(i-1,k)^2/a(i-1,k)-g*a(i-1,k)*h(i-1,k)/2)-.5*K(1,1)*a(i,k)*(q(i,k)^2/a(i,k)^2);
            up(i,k)=qp(i,k)/ap(i,k);
            hp(i,k)=ap(i,k)/b(i,1);
        end
            ap(1,k)=a(1,k);
            qp(1,k)=q(1,k);                                                     % boundary condition
            up(1,k)=u(1,k);
            hp(1,k)=h(1,k);
    
 
%%%%%%%   1st case d/s boundary condition by extrapolation
    
            ap(24,k)=ap(24-1,k);
            qp(24,k)=qp(24-1,k);
            up(24,k)=up(24-1,k);
            hp(24,k)=hp(24-1,k);
 
  
    for i=2:24-1
        ac(i,k)=a(i,k)-(delt/(delx))*(qp(i+1,k)-qp(i,k));
        qc(i,k)=q(i,k)+delt*g*ap(i,k)*(sx(i,1)-(n(1,1)^2*qp(i,k)^2*(b(i,1)+2*hp(i,k))^(4/3)*ap(i,k)^(-10/3)))-(delt/delx)*(qp(i+1,k)^2/ap(i+1,k)+g*ap(i+1,k)*hp(i+1,k)/2-qp(i,k)^2/ap(i,k)-g*ap(i,k)*hp(i,k)/2)-0.5*K(1,1)*ap(i,k)*(qp(i,k)^2/ap(i,k)^2);
        uc(i,k)=qc(i,k)/ac(i,k);
        hc(i,k)=ac(i,k)/b(i,1);
    end
        a(2:24-1,k+1)=(ap(2:24-1,k)+ac(2:24-1,k))/2;                        % boundary condition (B. C at u/s is provided at the begining)
        q(2:24-1,k+1)=(qp(2:24-1,k)+qc(2:24-1,k))/2;
        
            q(24,k)=0.4*q(1,k)+qp(24,k);                                    % intermediate boundary condition for tributaries. First tributary= 40% of inflow (RHEP).
                                                                            % This code was prepared for RHEP with two tributary contribution. For no tributary user can simply enter 0 in place of 0.4              
 
 %First Tributary          
           
        for i=25:62-1
            ap(i,k)=a(i,k)-(delt/(delx))*(q(i,k)-q(i-1,k));
            qp(i,k)=q(i,k)+delt*g*a(i,k)*(sx(i,1)-(n(1,1)^2*q(i,k)^2*(b(i,1)+2*h(i,k))^(4/3)*a(i,k)^(-10/3)))-(delt/delx)*(q(i,k)^2/a(i,k)+g*a(i,k)*h(i,k)/2-q(i-1,k)^2/a(i-1,k)-g*a(i-1,k)*h(i-1,k)/2)-.5*K(1,1)*a(i,k)*(q(i,k)^2/a(i,k)^2);
            up(i,k)=qp(i,k)/ap(i,k);
            hp(i,k)=ap(i,k)/b(i,1);
        end
            ap(1,k)=a(1,k);
            qp(1,k)=q(1,k);                                                 % boundary condition
            up(1,k)=u(1,k);
            hp(1,k)=h(1,k);
 
%%%%%%%   1st case d/s boundary condition by extrapolation
    
            ap(62,k)=ap(62-1,k);
            qp(62,k)=qp(62-1,k);
            up(62,k)=up(62-1,k);
            hp(62,k)=hp(62-1,k);
  
 
  
    for i=25:62-1
        ac(i,k)=a(i,k)-(delt/(delx))*(qp(i+1,k)-qp(i,k));
        qc(i,k)=q(i,k)+delt*g*ap(i,k)*(sx(i,1)-(n(1,1)^2*qp(i,k)^2*(b(i,1)+2*hp(i,k))^(4/3)*ap(i,k)^(-10/3)))-(delt/delx)*(qp(i+1,k)^2/ap(i+1,k)+g*ap(i+1,k)*hp(i+1,k)/2-qp(i,k)^2/ap(i,k)-g*ap(i,k)*hp(i,k)/2)-0.5*K(1,1)*ap(i,k)*(qp(i,k)^2/ap(i,k)^2);
        uc(i,k)=qc(i,k)/ac(i,k);
        hc(i,k)=ac(i,k)/b(i,1);
    end
         a(24:62-1,k+1)=(ap(24:62-1,k)+ac(24:62-1,k))/2;                     % boundary condition (B. C at u/s is provided at the begining)
        q(24:62-1,k+1)=(qp(24:62-1,k)+qc(24:62-1,k))/2;

           q(62,k)=0.3*q(1,k)+qp(62,k);                                     % 30% of inflow contributed by second tributary
                                     
  % Second Tributary                                                        % for no tributary simply put 0 in place of 0.3
    
         for i=63:YY-1
            ap(i,k)=a(i,k)-(delt/(delx))*(q(i,k)-q(i-1,k));
            qp(i,k)=q(i,k)+delt*g*a(i,k)*(sx(i,1)-(n(1,1)^2*q(i,k)^2*(b(i,1)+2*h(i,k))^(4/3)*a(i,k)^(-10/3)))-(delt/delx)*(q(i,k)^2/a(i,k)+g*a(i,k)*h(i,k)/2-q(i-1,k)^2/a(i-1,k)-g*a(i-1,k)*h(i-1,k)/2)-.5*K(1,1)*a(i,k)*(q(i,k)^2/a(i,k)^2);
            up(i,k)=qp(i,k)/ap(i,k);
            hp(i,k)=ap(i,k)/b(i,1);
         end
            ap(1,k)=a(1,k);
            qp(1,k)=q(1,k);                                                 % boundary condition
            up(1,k)=u(1,k);
            hp(1,k)=h(1,k);
 
%%%%%%%   1st case d/s boundary condition using sensor data
    
              if hp(YY-1,k)<h(S4,1)
                  hp(YY,k)=h(S4,1);
              else 
                  hp(YY,k)=hp(YY-1,k);
              end
        
              qp(YY,k)=qp(YY-1,k);  
              ap(YY,k)=hp(YY,k)*b(YY,1);

    for i=63:YY-1
        ac(i,k)=a(i,k)-(delt/(delx))*(qp(i+1,k)-qp(i,k));
        qc(i,k)=q(i,k)+delt*g*ap(i,k)*(sx(i,1)-(n(1,1)^2*qp(i,k)^2*(b(i,1)+2*hp(i,k))^(4/3)*ap(i,k)^(-10/3)))-(delt/delx)*(qp(i+1,k)^2/ap(i+1,k)+g*ap(i+1,k)*hp(i+1,k)/2-qp(i,k)^2/ap(i,k)-g*ap(i,k)*hp(i,k)/2)-0.5*K(1,1)*ap(i,k)*(qp(i,k)^2/ap(i,k)^2);
        uc(i,k)=qc(i,k)/ac(i,k);
        hc(i,k)=ac(i,k)/b(i,1);
    end
        a(62:YY-1,k+1)=(ap(62:YY-1,k)+ac(62:YY-1,k))/2;       % boundary condition (B. C at u/s is provided at the begining)
        q(62:YY-1,k+1)=(qp(62:YY-1,k)+qc(62:YY-1,k))/2;
        

    
        %%%%% Total Variation Diminishing 
    for i=2:YY-1
         rxp(i)=((a(i,k)-a(i-1,k))*(a(i+1,k)-a(i,k))+(a(i,k)*q(i,k)-a(i-1,k)*q(i-1,k))*(a(i+1,k)*q(i+1,k)-a(i,k)*q(i,k)))/...
                         ((a(i+1,k)-a(i,k))*(a(i+1,k)-a(i,k))+(a(i+1,k)*q(i+1,k)-a(i,k)*q(i,k))*(a(i+1,k)*q(i+1,k)-a(i,k)*q(i,k)));
         rxn(i)=((a(i,k)-a(i-1,k))*(a(i+1,k)-a(i,k))+(a(i,k)*q(i,k)-a(i-1,k)*q(i-1,k))*(a(i+1,k)*q(i+1,k)-a(i,k)*q(i,k)))/...
                          ((a(i,k)-a(i-1,k))*(a(i,k)-a(i-1,k))+(a(i,k)*q(i,k)-a(i-1,k)*q(i-1,k))*(a(i,k)*q(i,k)-a(i-1,k)*q(i-1,k)));           
    end
    
    for i=2:YY-1
       
                    rxp(i-1)=rxp(i);
                    rxn(i-1)=rxn(i);
                    rxp(1)=rxp(2);
                    rxn(1)=rxn(2);
    end
       
    for i=2:YY-1
                    rxp(i+1)=rxp(i);
                    rxn(i+1)=rxn(i);
                    rxp(YY)=rxp(YY-1);
                    rxn(YY)=rxn(YY-1);
    end
            
    for i=2:YY-1
                fxp(i)=max(0,min(2*rxp(i),1));
                fxn(i)=max(0,min(2*rxn(i),1));                      
    end
    
    for i=2:YY-1
      
                fxpn(i)=max(0,min(2*rxp(i-1),1));
                fxnp(i)=max(0,min(2*rxn(i+1),1));              
               
    end
        
    for i=2:YY-1
                Cr(i)=(u(i,k)+sqrt(g*h(i,k)))*delt/delx;
    end
          
            
    for i=2:YY-1
   
           if Cr(i)<=0.5
                    psi(i)=Cr(i)*(1-Cr(i));
           else
                    psi(i)=0.25;
           end
    end  
           
    for i=2:YY-1
                Gxp(i)=0.5*psi(i)*(1-fxp(i));
                Gxn(i)=0.5*psi(i)*(1-fxn(i));
                Gxpn(i)=0.5*psi(i)*(1-fxpn(i));
                Gxnp(i)=0.5*psi(i)*(1-fxnp(i));
                
    end
        
    for i=2:YY-1
             a(i,k+1)= a(i,k+1)      +(Gxp(i)+Gxnp(i))*(a(i+1,k)-a(i,k))-(Gxpn(i)+Gxn(i))*(a(i,k)-a(i-1,k)); 
             q(i,k+1)= q(i,k+1)      +(Gxp(i)+Gxnp(i))*(q(i+1,k)-q(i,k))-(Gxpn(i)+Gxn(i))*(q(i,k)-q(i-1,k));
                
    end   

     
        %%% 1st  case d/s boundary condition by extrapolation
%             
              
    
            a(YY,k+1)=a(YY-1,k+1); 
            q(YY,k+1)=q(YY-1,k+1);
            h(S4,k+1)=h(S4,k);

     % % correction for hydr0graph
% 
     for i=YY
            if q(i,k+1)<=q(1,1)
             q(i,k+1)=q(1,1);
             a(i,k+1)=a(i,k);
            end
     end 

        u(1:YY,k+1)=q(1:YY,k+1)./a(1:YY,k+1);
        h(1:YY,k+1)=a(1:YY,k+1)./b(1:YY,1);

        
       %%%% AREA VELOCITY IN NON CONSERVATIVE FORM
      
     for i=2:YY-1
            kvxh(i,k+1)=abs(h(i+1,k+1)-2*h(i,k+1)+h(i-1,k+1))/(abs(h(i+1,k+1))+abs(2*h(i,k+1))+abs(h(i-1,k+1)));
     end
     
     for i=1:YY
        kvxh(1,k+1)= abs(h(2,k+1)- h(1,k+1))/(abs(h(2,k+1))+abs(h(1,k+1)));
        kvxh(YY,k+1)=abs(h(YY,k+1)-h(YY-1,k+1))/(abs(h(YY,k+1))+abs(h(YY-1,k+1)));
     end
     
     for i=2:YY-1
            dvxhp(i,k+1)=(cons)*max(kvxh(i+1,k+1),kvxh(i,k+1));
            dvxhn(i,k+1)=(cons)*max(kvxh(i-1,k+1),kvxh(i,k+1));
     end
    
     for i=2:YY-1
            h(i,k+1)=h(i,k+1)+dvxhp(i,k+1)*(h(i+1,k+1)-h(i,k+1))-dvxhn(i,k+1)*(h(i,k+1)-h(i-1,k+1));
            u(i,k+1)=u(i,k+1)+dvxhp(i,k+1)*(h(i+1,k+1)-h(i,k+1))-dvxhn(i,k+1)*(h(i,k+1)-h(i-1,k+1));
     end
    
    
    for i=1:YY
        Surf(i,k+1)=h(i,k+1)+Ba(i);
    end
        Surf(YY,k+1)=Surf(YY-1,k+1);

      for i=YY-1:-1:1
            if Surf(i+1,k+1)>Surf(i,k+1)
               Surf(i,k+1)=Surf(i+1,k+1);
               hn(i,k+1)=Surf(i,k+1)-Ba(i,1);
               un(i,k+1)=(h(i,k+1)*u(i,k+1))/hn(i,k+1);
               h(i,k+1)=hn(i,k+1);
               u(i,k+1)=un(i,k+1);
               Surf(i,k+1)=h(i,k+1)+Ba(i,1); 
            end
       end
  
   
    for i=1:YY
        if h(i,k+1)>=hmax(i,1)
            hmax(i,1)=h(i,k+1);
            v_hmax(i,1)=u(i,k+1);
        end
    end
    
    for i=1:YY
        if u(i,k+1)>=vmax(i,1)                                             % maximum velocity
            vmax(i,1)=u(i,k+1);
    
        end
    end
      
    
    for i=1:YY                                                              % maximum surface
        Surfmax(i,k+1)=hmax(i,1)+Ba(i);
    end
    

   
      us(k)=q(1,k);                                                        % upstream hydrograph
      ds(k)=q(YY,k);                                                       % downstream hydrograph
      level(k)=Surf(YY,k);            
      NTq(k)=q(94,k);                                                      % Discharge in the downstream gauge section
      NTvel(k)=u(94,k);
      NTh(k)=h(94,k);
      NTlevel(k)=h(94,k)+91.2;                                             
      Depth_ds(k)=h(YY,k);
    
   
             
              
             CN= max(((u(1:YY,k)+sqrt(g*h(1:YY,k)))*delt/delx));            % Courant Number
              stst=max(abs(h(1:YY,k+1)-h(1:YY,k))./h(1:YY,k));              % Steady State
      
 end
 for t=1:z
   input(t)=us((t*3600)/delt)
   output(t)=ds((t*3600)/delt)
   Level(t)=level((t*3600)/delt)
   NT_vel(t)=NTvel((t*3600)/delt)                                          %Velocity at NT Road
   NT_depth(t)=NTh((t*3600)/delt);
   depth_p(t)= Depth_ds((t*3600)/delt);
   NT_Level(t)= NTlevel((t*3600)/delt);                                    % water level at NT Road Lakhimpur
   NT_Q(t)=NTq((t*3600)/delt)                                              %discharge at NT Road
 end
input=input';
NT_Q=NT_Q';
NT_vel=NT_vel';
NT_depth=NT_depth';
NT_Level=NT_Level';
end
