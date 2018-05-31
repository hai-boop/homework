%% @author Yvette-suyu
%% 2018.5.31
%% awesome ekf for free-fall


close all;
clear all;
Ts=1;                     %量测频率
t = 32;
N = fix(t/Ts);               %采样次数 
g = 9.81;
a = 100;
mrad=0.001;
X=zeros(2,N);                %目标真实位置、速度
X(:,1)=[5000;0];              %目标的初始位置，速度
Z=zeros(1,N);                %传感器对角度的测量
Q=diag([10.5,1]);              %过程噪声方差
R = 10;                       %量测噪声方差
F=[1 -Ts;0 1];                %状态转移矩阵F(k)
G=[-0.5*Ts.^2;1];             %输入矩阵G(k)
h=@(x)(atan(x/a));                      % 观测矩阵 H(k) theta与速度无关 转化为一维向量
W = sqrt(Q)*randn(2,1);      % 模拟产生过程噪声(高斯分布的随机噪声)
V = sqrt(R)*mrad*randn(1,1);           % 模拟产生测量噪声
 
%假设 物体下落时与原点的水平距离为100m

%% real trace
for k=2:N
    X(:,k)=F*X(:,k-1)+G*g+W;%实际的位置（高度）与速度
end

for k=1:N

    Z(k)=atan(X(1,k)/a)+V;
%     Z(k)=h(X(1,k))+V;%实际的角度，且角度变化由高度表示
end

%% EKF
Xekf=zeros(2,N);
Xekf(:,1)=X(:,1);  %EKF初始化
P0=eye(2);          %误差协方差矩阵初始化
z=zeros(1,N);
Ph=zeros(1,N);
Ph(1)=P0(1,1);
Pv=zeros(1,N);
Pv(1)=P0(2,2);
Pm=zeros(1,N);

for k=2:N
    Xn = F*Xekf(:,k-1)+G*g; %一步预测
    P1 = F*P0*F'+Q;     %预测误差协方差
   % ob = atan(Xn/a);    %观测预测
    z(k-1) = h(Xn(1)); 
   % 求解雅可比矩阵H
    Hh = JacobianH(Xn(1));
   
%     S = Hh*P1*Hh'+R;
%     K = P1*Hh'*S^-1;
    K = P1*Hh'*inv(Hh*P1*Hh'+R);%EKF增益
    Xekf(:,k) = Xn+K*(Z(k-1)-z(k-1));%状态更新
    P0 = (1-K*Hh)*P1;%滤波误差协方差更新
    Ph(k)=P1(1,1);
    
end

%%误差分析
Err_Messure=zeros(1,N);
Err_EKFh=zeros(1,N);
Err_EKFv=zeros(1,N);
for k=1:N
    Err_Messure(k)=abs(cos(z(k))-cos(Z(k)));

    Err_EKFv(k)=sqrt(0.5*(Xekf(2,k)-X(2,k))^2);
end
%% figure

t=1:N;
figure
plot(t,X(1,:),'-b')
hold on;
plot(t,Xekf(1,:),'-ks');
legend('real','ekf extimate');         
xlabel('time/s');
ylabel('H/m');
title('EKF Simulation of Height');

figure
plot(t,X(2,:),'-b');
hold on;
plot(t,Xekf(2,:),'-ks');
legend('real','ekf extimate');         
xlabel('time/s');
ylabel('v/m/s');
title('EKF Simulation of Velocity');

figure
plot(t,cos(Z),'-b');
hold on;
plot(t,cos(z),'-ks');
legend('real','ekf extimate');         
xlabel('time/s');
title('EKF Simulation of cos angle');

figure
plot(t,Err_Messure,'-b');
legend('messure error');         
xlabel('time/s');

figure
plot(t,Ph,'-b');
legend('Height error');         
xlabel('time/s');

figure
plot(t,Err_EKFv,'-b');
legend('Velocity error');         
xlabel('time/s');

