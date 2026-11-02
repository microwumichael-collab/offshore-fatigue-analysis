clear
clc
close all

g = 9.81;
d = 63.2;
dSWL = 3.16;
d = d+dSWL;
D_z0 = 8.18;
rho = 1025;
D = 11.08;
A = 1/4*pi*D^2;
Us = 1.56;
H = 6.7;

T = 14.3; 


omega = 2*pi/T;

k = fDispersionV5(d,T,H,2); 
Ca = 1;
Cd = 9.52E-01;

%% All calculation done assuming x=t=0
eta_max = H/2+H^2*k/16 *cosh(k*d)/sinh(k*d)^3*(2+cosh(2*k*d));
z = linspace(0,d+eta_max,100);
dz = z(2)-z(1);

u = H*9.81*k*cosh(k*(z))/(2*omega*cosh(k*d))+3/16*H^2*omega*k*cosh(2*k*(z))/sinh(k*d)^4;
ut = @(t) H*9.81*k*cosh(k*(z))/(2*omega*cosh(k*d))*cos(-omega*t)+3/16*H^2*omega*k*cosh(2*k*(z))/sinh(k*d)^4*cos(-2*omega*t);

U = H*9.81*k*cosh(k*(d))/(2*omega*cosh(k*d))+3/16*H^2*omega*k*cosh(2*k*(d))/sinh(k*d)^4;
Uc = Us*(z/d).^(1/7);
a = cosh(k*z)/cosh(k*d);
b = cosh(2*k*z)/sinh(k*d)^4;
c = (z/d).^(1/7);

E1 = sum(a.^2)*dz;
E2 = sum(a.*b)*dz;
E3 = sum(a)*dz;
E4 = sum(b.^2)*dz;
E5 = sum(b)*dz;
Eu1 = sum(a.*c)*dz;
Eu2 = sum(b.*c)*dz;
Eu3 = sum(c.^2)*dz;

F = @(t) 0.5*rho*Cd*D*(H^2/4*g^2*k^2/omega^2*E1*cos(-omega*t)^2+3/16*H^3*g*k^2*E2*cos(-omega*t)*cos(-2*omega*t)+H*g*k/omega*Eu1*cos(-omega*t)*Us+9/256*H^4*omega^2*k^2*E4*cos(-2*omega*t)^2+3/8*H^2*omega*k*Eu2*cos(-2*omega*t)*Us+Us^2*Eu3)+rho*(1+Ca)*A*(H/2*g*k*sin(-omega*t)*E3+3/8*H^2*omega^2*k*sin(2*(-omega*t))*E5);
syms t
F_expr = F(t);
dfdt = diff(F_expr, t);
dF = matlabFunction(dfdt, 'Vars', t);


tt = linspace(0, T, 1001);
s  = sign(dF(tt));
ix = find(s(1:end-1) ~= s(2:end));   
cp = arrayfun(@(i) fzero(dF, [tt(i) tt(i+1)]), ix);
F_at_cp = arrayfun(F, cp);

[Fmax, i] = max(F_at_cp);
t_max = cp(i);

[Fmin, i_min] = min(F_at_cp);
t_min = cp(i_min);

%% Overturning Moment
h = 0:dz:63.2+19;
u_tmax = H*9.81*k*cosh(k*(h))/(2*omega*cosh(k*d))*cos(-omega*t_max)+3/16*H^2*omega*k*cosh(2*k*(h))/sinh(k*d)^4*cos(-2*omega*t_max)+Us*(h/d).^(1/7);
v_tmax = H*9.81*k*cosh(k*(h))/(2*cosh(k*d))*sin(-omega*t_max)+3/8*H^2*omega^2*k*cosh(2*k*(h))/sinh(k*d)^4*sin(-2*omega*t_max);
f = rho*(1+Ca)*A*v_tmax+0.5*rho*Cd*D*u_tmax.*abs(u_tmax);
f(h>max(z)) = 0;
f(end) = 31e5/dz;

BMD = zeros(1, length(f));
V = zeros(1, length(f));

for i2=1:length(f)
    BMD(i2) = sum(f(i2:length(f)).*(h(i2:length(h))-h(i2)))*dz+357e6;
    V(i2) = sum(f(i2:end))*dz;
end
%
FHmax_val = F(t_max);
Vmax_val = max(V);
BMDmax_val = max(BMD);

u_tmin = H*9.81*k*cosh(k*(h))/(2*omega*cosh(k*d))*cos(-omega*t_min)+3/16*H^2*omega*k*cosh(2*k*(h))/sinh(k*d)^4*cos(-2*omega*t_min)+Us*(h/d).^(1/7);
v_tmin = H*9.81*k*cosh(k*(h))/(2*cosh(k*d))*sin(-omega*t_min)+3/8*H^2*omega^2*k*cosh(2*k*(h))/sinh(k*d)^4*sin(-2*omega*t_min);


eta_tmin = H/2 * cos(-omega*t_min) + H^2*k/16 * cosh(k*d)/sinh(k*d)^3 * (2+cosh(2*k*d)) * cos(-2*omega*t_min);
current_water_level = d + eta_tmin; 

f_min = rho*(1+Ca)*A*v_tmin + 0.5*rho*Cd*D*u_tmin.*abs(u_tmin);


f_min(h > current_water_level) = 0; 
f_min(end) = 31e5/dz; 

BMD_min = zeros(1, length(f_min));
V_min = zeros(1, length(f_min));


for i2=1:length(f)

    BMD(i2) = sum(f(i2:length(f)).*(h(i2:length(h))-h(i2)))*dz + 357e6;
    V(i2) = sum(f(i2:end))*dz;

    BMD_min(i2) = sum(f_min(i2:length(f_min)).*(h(i2:length(h))-h(i2)))*dz + 357e6; 
    V_min(i2) = sum(f_min(i2:end))*dz;
end






delta_shear = V - V_min; 
delta_BMD = BMD - BMD_min;


h1 = [4.76
2.5
-0.5
-4
-7.5
-11.06
-14.62
-18.18
-21.74
-25.3
-28.86
-32.381
-36.261
-40.141
-44.021
-47.901
-51.781
-55.661
-59.541
-63.2
-67.08
-70.96
-74.84
-78.72
-82.6
-84.6
-86.6
-89.1
0
18.7
15.5
12
8.5
5.3
4.85
1
-3]+63.2;


Delta_BMD_at_h1 = interp1(h, delta_BMD, h1, 'spline');
Delta_shear_at_h1=interp1(h, delta_shear, h1, 'spline');
Delta_BMD_at_h1(Delta_BMD_at_h1<0) = 0; 
Delta_shear_at_h1(Delta_shear_at_h1<0) = 0;
%% --- 绘制 Delta 分布图 ---
figure
subplot(1,2,1)
hold on
% 绘制剪力幅值图
plot(delta_shear, h, 'LineWidth', 2, 'Color', '#D95319'); % 使用显眼的橙红色
title('Shear Force Range (\Delta V)')
ylabel('depth [m]')
xlabel('\Delta V [N]')
% 根据实际数值范围调整 xlim，这里先让它自适应或给一个宽泛的范围
% xlim([0, max(delta_shear)*1.1]) 
grid on
hold off

subplot(1,2,2);
hold on
% 绘制弯矩幅值图
plot(delta_BMD, h, 'LineWidth', 2, 'Color', '#7E2F8E'); % 使用紫色
title('Bending Moment Range (\Delta BMD)')
ylabel('depth [m]')
xlabel('\Delta BMD [N\cdot m]')
% xlim([0, max(delta_BMD)*1.1])
grid on
hold off

sgtitle(['Fatigue Load Ranges (Fixed T = ', num2str(T), ' s)'])
%% --- 绘制剪力与弯矩包络图 (Envelope Diagram) ---
figure
% 1. 剪力对比图
subplot(1,2,1)
hold on
plot(V, h, 'LineWidth', 2, 'Color', 'b');       % 最大剪力 (波峰)
plot(V_min, h, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--'); % 最小/反向剪力 (波谷)
xline(0, 'k-', 'LineWidth', 1);                 % 添加 0 刻度基准线
title('Shear Force Envelope')
ylabel('Depth [m]')
xlabel('V [N]')
legend('Maximum (Wave Crest)', 'Minimum (Wave Trough)', 'Location', 'best')
grid on
hold off

% 2. 弯矩对比图
subplot(1,2,2);
hold on
plot(BMD, h, 'LineWidth', 2, 'Color', 'b');     % 最大弯矩 (波峰)
plot(BMD_min, h, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--'); % 最小弯矩 (波谷)
xline(0, 'k-', 'LineWidth', 1);                 % 添加 0 刻度基准线
title('Bending Moment Envelope')
ylabel('Depth [m]')
xlabel('Bending Moment [N\cdot m]')
legend('Maximum (Wave Crest)', 'Minimum (Wave Trough)', 'Location', 'best')
grid on
hold off

sgtitle(['Load Envelopes for Fixed T = ', num2str(T), ' s'])