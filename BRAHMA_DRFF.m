z=23;
A_S=zeros([z,1]);
A_S1=zeros([z,1]);
power_met=zeros([z,1]);

inflow= xlsread('Inflow.xlsx',1);
n2=24;
n1=zeros([z,1]);
n1(5:(5+n2-1),1)=1;
n1(29:(29+n2-1),1)=1;
n1(53:(53+n2-1),1)=1;
safe_discharge=250*3600*10^-6;
advanced_spill=zeros([z 1]);
lead_time=input('Enter lead time: ');                                       % we used 6 hrs as per the position of the upstream gauged site
lead_time=6;
WL=input('Enter Warning Level at d/s: ');                                   % we used 93.8 m as WL for the dowstream flood prone section
st=60;
FM=input('Enter minimum power: ');                                          % the power for 1 unit in this case is 135 MW
PC=input('Enter Plant Capacity: ');                                         % In this case the plant capacity is 405 MW
beta=input('Enter beta: ');                                                 % we used the pre-release factor beta as 0.61

[f_s,NT_Level,p_p,a_s,r_r,NT_Q,NT_h,NT_vel,CN,output]=BRAHMA_DRFF_func(inflow,advanced_spill,z,st,n1)
[Qmax,idx]=max(inflow);
for s=1:z
    if NT_Level(s,1)>WL
          if NT_Level(s-1,1)<WL
             index=idx;
             back_index=index-lead_time;
             for r =back_index:index-1
                 qi(r,1)=NT_Q(r,1);
                 A_S(r,1)=beta*safe_discharge;
                 
             end
          end
    end
end

[f_s1,NT_Level1,p_p1,a_s1,r_r1,NT_Q1,NT_h1,NT_vel1,CN1,output1]= BRAHMA_DRFF_func(inflow,A_S,z,st,n1);

% Generate performance metrics arrays
FR_SOP = zeros(z, 1);       
PR_SOP= zeros(z, 1);
FR_AOP =zeros(z, 1);       
PR_AOP= zeros(z, 1);

for u_i = 1:z
   % Determine flood reliability for SOP
    if NT_Level(u_i) > WL
        FR_SOP(u_i) = 0;
    else
        FR_SOP(u_i) = 1;
    end

    % Determine flood reliability for AOP
    if NT_Level1(u_i) > WL
        FR_AOP(u_i) = 0;
    else
        FR_AOP(u_i) = 1;
    end

  % Determine power met for AOP
   if p_p(u_i) < FM 
       PR_SOP(u_i) = 0;
   elseif p_p(u_i) == FM 
       PR_SOP(u_i) = 0.33;
   elseif p_p(u_i) == 2 * FM 
       PR_SOP(u_i) = 0.66;
   elseif p_p(u_i) == PC
       PR_SOP(u_i) = 1;
   end

    % Determine power met for AOP
   if p_p1(u_i) < FM 
       PR_AOP(u_i) = 0;
   elseif p_p1(u_i) == FM 
       PR_AOP(u_i) = 0.33;
   elseif p_p1(u_i) == 2 * FM 
       PR_AOP(u_i) = 0.66;
   elseif p_p1(u_i) == PC 
       PR_AOP(u_i) = 1;
   end
end

% Compute flood count
FC_SOP = z-(sum(FR_SOP));
FC_AOP= z-(sum(FR_AOP));

Improved_FC=FC_SOP-FC_AOP;

% Concatenate values
values = [r_r(:,1), p_p(:,1), a_s(:,1), NT_Q(:,1), ...
          NT_h(:,1), NT_Level(:,1), NT_vel(:,1), ...
          PR_SOP, FR_SOP]; 

values_adj = [r_r1(:,1), p_p1(:,1), a_s1(:,1), NT_Q1(:,1), ...
              NT_h1(:,1), NT_Level1(:,1), NT_vel1(:,1), ...
              PR_AOP, FR_AOP];

headers = {'Release', 'Power Generated', 'Active Storage', 'Downstream Discharge', 'Downstream Depth', 'Downstream Level', 'Downstream Velocity', 'Power Met', 'Reliability'};


% Combine headers with values
values_with_headers = [headers; num2cell(values)];       % Add headers to SOP results
values_with_headers1 = [headers; num2cell(values_adj)];  % Add headers to AOP results

writecell(values_with_headers, 'Result_SOP.xlsx'); 
writecell(values_with_headers1, 'Result_AOP.xlsx');


xlswrite('Result_SOP.xlsx', values_with_headers); 
xlswrite('Result_AOP.xlsx', values_with_headers1); 

 plot(1:z,NT_Level,1:z,NT_Level1)
 title('Levels')
 xlabel 'Time(hrs)';
 ylabel 'Discharge (hrs)';
 legend('NT_Level','NT_Level1');
