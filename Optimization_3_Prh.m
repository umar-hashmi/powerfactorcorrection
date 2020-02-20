clear
close all
clc
load('real_time_data')
tic
pts = 96; %%number of points considered
time=real(1:pts,2);

global P_g_i Q_g_i S_B_max

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
b_int=1000;
B_0 = b_int*ones(pts,1);

h=0.25;

S_B_max = 0.9*del_max*h/e_ch;



for i=1:96
    price=real(i:pts,1);
    C_cap =S_B_max^2;
    x_upper= del_max*h*ones(pts-i+1,1);
    x_lower= del_min*h*ones(pts-i+1,1);
    B_max = 2000*ones(pts - i+1,1);
    B_min = 200*ones(pts-i+1,1);
    M = tril(ones(pts-i+1,pts-i+1));
    if i >1
        B_0 = t1*ones(pts-i+1,1);
    end
    tx=0;
    P_g_i =P_g(i);

    Q_g_i =Q_g(i);
    i
    
    cvx_begin
    variables x_ch(pts-i+1,1) x_ds(pts-i+1,1) b_q(pts-i+1,1)
    variables b_p(pts-i+1,1) b_p_bat(pts-i+1,1)
    minimize sum(price'*x_ch/e_ch-price'*x_ds*e_dis)     %trace((x_ch/e_ch-x_ds*e_dis)*price')
    subject to
    zeros(pts-i+1,1)<= x_ch <= x_upper;    %%charging
    zeros(pts-i+1,1)<= x_ds <= -x_lower;    %% discharging
    b_p_bat == x_ch-x_ds;               %%change in charge level of battery
    b_p == x_ch/e_ch-x_ds*e_dis;        %%active power output of battery
    B_min <= B_0 + M*b_p_bat <= B_max;  %%battery capacity
    (b_p).^2 + (b_q).^2 <= C_cap;   %%apparent power limit
    - tan_limit.*(b_p(1)+P_g(i)) + (b_q(1)+Q_g(i)) <= 0;  %%power factor
    - tan_limit.*(b_p(1)+P_g(i)) - (b_q(1)+Q_g(i)) <= 0;  %%power factor
    b_p(1)+P_g(i)>=0
    cvx_end
    
    temp1 =sum(price'*x_ch/e_ch-price'*x_ds*e_dis) ; P_temp1 =b_p(1); Q_temp1 = b_q(1);
    t_act1 = isnan(b_p_bat(1));
    
    cvx_begin
    variables x_ch(pts-i+1,1) x_ds(pts-i+1,1) b_q(pts-i+1,1)
    variables b_p(pts-i+1,1) b_p_bat(pts-i+1,1)
    minimize sum(price'*x_ch/e_ch-price'*x_ds*e_dis)     %trace((x_ch/e_ch-x_ds*e_dis)*price')
    subject to
    zeros(pts-i+1,1)<= x_ch <= x_upper;    %%charging
    zeros(pts-i+1,1)<= x_ds <= -x_lower;    %% discharging
    b_p_bat == x_ch-x_ds;               %%change in charge level of battery
    b_p == x_ch/e_ch-x_ds*e_dis;        %%active power output of battery
    B_min <= B_0 + M*b_p_bat <= B_max;  %%battery capacity
    (b_p).^2 + (b_q).^2 <= C_cap;   %%apparent power limit
    - tan_limit.*(b_p(1)+P_g(i)) + (b_q(1)+Q_g(i)) >= 0;  %%power factor
    - tan_limit.*(b_p(1)+P_g(i)) - (b_q(1)+Q_g(i)) >= 0;  %%power factor
    b_p(1)+P_g(i)<=0
    cvx_end
     
    temp2 =sum(price'*x_ch/e_ch-price'*x_ds*e_dis) ; P_temp2 =b_p(1); Q_temp2 = b_q(1);
    t_act2 = isnan(b_p_bat(1));
    
    if t_act1 == 0 && t_act2 == 0
        if temp1> temp2
            P_B(i) = P_temp2;
            Q_B(i) = Q_temp2;
        else
            P_B(i) = P_temp1;
            Q_B(i) = Q_temp1;
        end
    elseif t_act1 == 0 && t_act2 == 1
        P_B(i) = P_temp1;
        Q_B(i) = Q_temp1;
    elseif t_act1 == 1 && t_act2 == 0
        P_B(i) = P_temp2;
        Q_B(i) = Q_temp2;
    else
        tan_temp = Q_g(i)/P_g(i);
        if tan_temp >= tan_limit && P_g(i) >= 0
            lb_lim = 0;
            ub_lim = min((B_max(1) - bcap(i-1))/e_ch, del_max*h/e_ch);
        elseif tan_temp >= tan_limit && P_g(i) < 0
            lb_lim = max((B_min(1) - bcap(i-1))*e_dis, del_min*h*e_dis);
            ub_lim = 0;
        elseif tan_temp <= -tan_limit && P_g(i) >= 0
            lb_lim = 0;
            ub_lim = min((B_max(1) - bcap(i-1))/e_ch, del_max*h/e_ch);
        else 
            lb_lim = max((B_min(1) - bcap(i-1))*e_dis, del_min*h*e_dis);
            ub_lim = 0;
        end
        
        P_B(i) = fminbnd(@tanreal,lb_lim ,ub_lim );
%         Q_B(i) = -sign(Q_g(i))*sqrt(S_B_max^2 - (P_B(i))^2);
        if Q_g(i) >= 0
            Q_B(i) = max(-sqrt(S_B_max^2 - (P_B(i))^2), -Q_g(i));
        else
            Q_B(i) = min(sqrt(S_B_max^2 - (P_B(i))^2), -Q_g(i));
        end
            
%         if P_g(i) >= 0
%             P_B(i) = min((B_max(1) - bcap(i-1))/e_ch, del_max*h/e_ch);
%             if abs(sqrt(S_B_max^2 - (P_B(i))^2)) > abs(Q_g(i))
%                 Q_B(i) = -Q_g(i);
%             else
%                 Q_B(i) = -sign(Q_g(i))*sqrt(S_B_max^2 - (P_B(i))^2);
%             end
%         else
%             P_B(i) = max((B_min(1) - bcap(i-1))*e_dis, del_min*h/e_dis);
%             if abs(sqrt(S_B_max^2 - (P_B(i))^2)) > abs(Q_g(i))
%                 Q_B(i) = -Q_g(i);
%             else
%                 Q_B(i) = -sign(Q_g(i))*sqrt(S_B_max^2 - (P_B(i))^2);
%             end
%         end           
    end
    x_charge = subplus(P_B(i))*e_ch - subplus(-P_B(i))/e_ch;
    
    t1=B_0(1,1) + x_charge;
    bcap(i,1) = B_0(1,1) + x_charge;
    clear B_0
    clear price
    clear x_lower
    clear x_upper
    clear x_ds
    clear x_ch
    clear b_q

end

xopt=bcap-[b_int; bcap(1:end-1)];
pr=real(1:pts,1);
profit=  sum(pr'*P_B')/1000

pf_cvx = abs(P_B'+P_g)./sqrt((P_B'+P_g).^2 + (Q_B'+Q_g).^2 );

figure; plot(time,pf_cvx)

a8=sort(pf_cvx);

mean(pf_cvx)

C_cap =S_B_max^2*ones(pts,1);
converter_usage_fac=mean(((P_B).^2 + (Q_B).^2)./S_B_max^2)


MINPF=min(pf_cvx)
