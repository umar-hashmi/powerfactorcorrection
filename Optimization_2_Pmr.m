clear 
close all
clc
load('real_time_data')
tic
pts = 96; %%number of points considered
time=real(1:pts,2);
price=real(1:pts,1);

% Suppose we set the lower bound of power factor
pf_limit=0.9;
theta = acos(pf_limit);
tan_limit = tan(theta);

% B_max=1; %Capacity of battery is 1 kWh
load('15_min_load_data.mat')
load('15_min_pv_data.mat')

P_g = avg_active_grid(1:pts,1);
Q_g = avg_reactive_grid(1:pts,1);

e_ch=0.95;
e_dis =0.95;
del_max = 4000;
del_min = -del_max;
B_0 = 1000*ones(pts,1);
B_max = 2000*ones(pts,1);
B_min = 200*ones(pts,1);
h=0.25;
x_upper= del_max*h*ones(pts,1);
x_lower= del_min*h*ones(pts,1);
S_B_max = 1.5*del_max*h/e_ch;
C_cap =S_B_max^2*ones(pts,1);
M = tril(ones(pts,pts));
 
cvx_begin 
cvx_precision best
cvx_solver Gurobi 
cvx_solver_settings( 'MIPGap', 1e-60 ) 
% cvx_solver_settings( 'Method', 1 ) 
cvx_solver_settings('MIPFocus',2) 
% cvx_solver Mosek 
% MIPFocus=2;
% cvx_solver_settings('eps', 1e-8)
variables x_ch(pts,1) x_ds(pts,1) b_q(pts,1) 
variables b_p(pts,1) b_p_bat(pts,1) y(pts,1)
variable z(pts,1) binary
minimize sum(price'*b_p)    %trace((x_ch/e_ch-x_ds*e_dis)*price')
subject to 
    zeros(pts,1)<= x_ch <= x_upper;    %%charging 
    zeros(pts,1)<= x_ds <= -x_lower;    %% discharging 
    b_p_bat == x_ch-x_ds;               %%change in charge level of battery
    b_p == x_ch/e_ch-x_ds*e_dis;        %%active power output of battery
    B_min <= B_0 + M*b_p_bat <= B_max;  %%battery capacity
    (b_p).^2 + (b_q).^2 <= C_cap;   %%apparent power limit
%     -2*tan_limit.*y + tan_limit.*(b_p+P_g) <= (b_q+Q_g) <= 2*tan_limit.*y - tan_limit.*(b_p+P_g); 
    + tan_limit.*(b_p+P_g) -2*tan_limit*y + (b_q+Q_g) <= zeros(pts,1);  %%power factor
    + tan_limit.*(b_p+P_g) -2*tan_limit*y - (b_q+Q_g) <= zeros(pts,1); %%power factor
    2*y - (b_p+P_g) >= 0;
    y >= z.*(P_g + x_lower);
    y >= (b_p+P_g) + z.*(P_g + x_upper/e_ch) - (P_g + x_upper/e_ch);
    y <= z.*(P_g + x_upper/e_ch);
    y <= (b_p+P_g) + z.*(P_g + x_lower) - (P_g + x_lower);
%     y >= (P_g +(B_min - B_0 - M*b_p_bat)*e_dis);
%     y <= (P_g +(B_max - B_0 - M*b_p_bat)/e_ch);
cvx_end
toc
profit=  sum(price'*b_p)/1000

B = B_0 + M*b_p_bat;

time_s = h: h: length(B)/4;

pf_cvx = abs(b_p+P_g)./sqrt((b_p+P_g).^2 + (b_q+Q_g).^2 );

figure; plot(time_s,pf_cvx)

x_opt = x_ch - x_ds;
pr=real(1:pts,1);

a8=sort(pf_cvx);

mean(pf_cvx)

MINPF=min(pf_cvx)

figure; plot(((b_p).^2 + (b_q).^2)./C_cap)

converter_usage_fac=sum(((b_p).^2 + (b_q).^2)./C_cap)/96

figure; plot(x_opt)
