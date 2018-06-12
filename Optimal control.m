%% Aircraft Longitudinal Motions - STOL Aircraft
clear; close all; clc;

%% 系统状态空间模型
A = [-0.0397 -0.280 0 -0.562;
    0.135 -0.538 -0.957 0;
    0.0207 0.441 -1.410 0
    0 0 1 0];

B_u = [-0.0052 -0.102;
    0.031 0.037;
    -1.46 -0.066;
    0 0];

B_w = [0.0397 0.280;
    -0.135 0.538;
    0.0207 -0.441;
    0 0];

C = [-0.0178 -1.92 0 1.92;
    1 0 0 0];

%% 系统初始状态
x(:,1) = [370; 14; 1; 14];
y(:,1) = C * x(:,1);
t(:,1) = 0;

%% 调参
v = 300;
u_w = 10;
alpha_w = u_w / v;
Q = [100 0 0 0;
   0 1 0 0;
   0 0 1 0;
   0 0 0 2.25];
R = [100 0;
   0 1];
delta = 0.1; % 采样时间
t_suml = 500; % 仿真时间

%% LQR问题求解
B = B_u * inv(R) * B_u';
P = are(A,B,Q); % 求解Riccati微分方程
K = inv(R) * B_u' * P; % Kalman增益
disp("Kalman增益： K = ");
disp(K);

%% 系统仿真迭代更新

% ---------- 施加最优控制 ----------
for k = 1 : t_suml/delta
    u_oc(:,k) = -K * x(:,k); % 最优控制量
    x1(:,k) = A * x(:,k) + B_u * u_oc(:,k) + B_w * [u_w; alpha_w];
    % x1(:,k) = A * x(:,k) + B_u * u_oc(:,k);
    x(:,k+1) = x(:,k) + x1(:,k) * delta;
    y_oc(:,k+1) = C * x(:,k+1);  
    t(:,k+1) = k * delta; % 采样时刻
end
% ---------- 未施加最优控制的控制量 ----------
% u = [-13.1581 * ones(1,t_suml/delta); 32.9386 * ones(1,t_suml/delta)]; % 控制量
u = [-30.9203 * ones(1,t_suml/delta); 421.2383 * ones(1,t_suml/delta)]; % 控制量
% ---------- 未施加最优控制（有扰动） ----------
for k = 1 : t_suml/delta
    x1(:,k) = A * x(:,k) + B_u * u(:,k) + B_w * [u_w; alpha_w];
    % x1(:,k) = A * x(:,k) + B_u * u(:,k);
    x(:,k+1) = x(:,k) + x1(:,k) * delta;
    y(:,k+1) = C * x(:,k+1);  
end

% ---------- 未施加最优控制（无扰动） ----------
for k = 1 : t_suml/delta
    x1(:,k) = A * x(:,k) + B_u * u(:,k);
    x(:,k+1) = x(:,k) + x1(:,k) * delta;
    y_without_disturbance(:,k+1) = C * x(:,k+1);  
end

%% 绘图
% ---------- Velocity on Altitude Change with Time ----------
figure(1);
plot(t, y_oc(1,:),'r');
hold on;
plot(t, y(1,:),'b');
hold on;
plot(t, y_without_disturbance(1,:),'g')
legend('最优控制','系统状态','未受扰动系统状态');
xlabel('time  (s)');
ylabel('Velocity on Altitude  (m/s)');
title('Velocity on Altitude Change with Time');

% ------ Velocity on longitudinal Displacement Change with Time -------
figure(2);
plot(t, y_oc(2,:),'r'); % 水平速度
hold on;
plot(t, y(2,:),'b');
hold on;
plot(t, y_without_disturbance(2,:),'g');
legend('最优控制','系统状态','未受扰动系统状态');
xlabel('time (s)');
ylabel('Velocity on longitudinal Displacement  (m/s)');
title('Velocity on longitudinal Displacement Change with Time');

% ---------- Optimal Control Rate Change with Time ----------
figure(3);
plot(t(1:t_suml/delta), u_oc(1,:),'r');
hold on;
plot(t(1:t_suml/delta), u_oc(2,:),'g');
legend('控制量1','控制量2');
xlabel('time (s)');
ylabel('Optimal Control Rate');
title('Optimal Control Rate Change with Time');

%% 计算性能指标
% 竖直方向
% ---------- 计算超调量 ----------
overshoot_h = (max(y(1,:)) - y(1,k)) / y(1,k); % 超调量
% ---------- 计算调整时间 ----------
for i = k : -1 : 1
    if abs(y(1,i) - y(1,k)) > 0.02 * y(1,k);
        break;
    end
end
ts_h = t(i-1); % 调整时间
% ---------- 计算延迟时间 ----------
for i = 1 : k
    if y(1,i) <= y(1,k) / 2 && y(1,i+1) >= y(1,k) / 2
        break;
    end
end
if i < 2000
    have_td_h = 1; % 标志位：等于1，表示有此项性能指标；等于0，表示无此项性能指标
    td_h = t(i); % 延迟时间
else
    have_td_h = 0;
end
% ---------- 显示性能指标 ----------
disp("竖直方向性能指标：");
fprintf('    超调量：%f%%\n',overshoot_h * 100);
fprintf('    调整时间：%fs\n',ts_h);
if have_td_h == 1
    fprintf('    延迟时间：%fs\n',td_h);
else if have_td_h == 0
        fprintf('    无延迟时间此项性能指标');
    end
end
% 纵向速度
% ---------- 计算超调量 ----------
overshoot_x = (max(y(2,:)) - y(2,k-1)) / y(2,k-1); % 超调量
% ---------- 计算调整时间 ----------
for j = k : -1 : 1
    if abs(y(2,j) - y(2,k)) > 0.02 * y(2,k);
        break;
    end
end
ts_x = t(j-1); % 调整时间
% ---------- 计算延迟时间 ----------
for j = 1 : k
    if y(2,j) <= y(2,k) / 2 && y(2,j+1) >= y(2,k) / 2
        break;
    end
end
if j < 2000
    have_td_x = 1;
    td_x = t(j); % 延迟时间
else
    have_td_x = 0;
end
% ---------- 显示性能指标 ----------
disp("纵向性能指标：");
fprintf('    超调量：%f%%\n',overshoot_x * 100);
fprintf('    调整时间：%fs\n',ts_x);
if have_td_x == 1
    fprintf('    延迟时间：%fs\n',td_x);
else if have_td_x == 0
        fprintf('    无延迟时间这项性能指标\n');
     end
end
