%% author: yvette-suyu
%% 2018.5.20

% Kalman
t=0.01;
N=32;
A=[1 -t;0 1];                               % 状态转移矩阵 Φ(k)
B=[-0.5*t.^2;1];
H=[1 0];                                    % 观测矩阵 H(k)
g=9.81;

X=zeros(2,N);
X(:,1)=[50000;0];                           % 目标的状态向量 X(k)
Z=zeros(1,N);
Z(1)=H*X(:,1);                              %观测量Z(k)  
P=zeros(N,2);
P0=[15 0;0 1];                            % 一步预测的协方差 P(k)
P(1,1)=P0(1,1);
P(1,2)=P0(2,2);

Xkf=zeros(2,N);
Xkf(:,1)=X(:,1);

Q = [0,0;0,0];
R =1;
W = sqrt(Q)*randn(2,N);                               % 模拟产生过程噪声(高斯分布的随机噪声)
V = sqrt(R)*randn(1,N);                               % 模拟产生测量噪声


for k=2:N

    X(:,k) = A * X(:,k-1)+g*B+V(k-1);               % 状态方程:X(k+1)=Φ(k)X(k)+G(k)V(k),其中G(k)=1
    Z(:,k)=H*X(:,k)+V(k);                           % 观测方程:Z(k+1)=H(k+1)X(k+1)+W(k+1),Z(k+1)是k+1时刻的观测值
end

% Q=std(V)^2;                                 % W(k)的协方差,std()函数用于计算标准偏差  
% R=std(W)^2;                                 % V(k)的协方差 covariance
                            

for t=2:N
    
   
    X_pred=A*Xkf(:,t-1)+g*B;
    P_pred=A*P0*A'+Q;                   % 一步预测的协方差 P(k+1|k)
    K=P_pred*H'*inv(H*P_pred*H'+R);     % 卡尔曼滤波器增益 K
    Xkf(:,t)=X_pred+K*(Z(t)-H*X_pred);  % 状态更新方程 X(k+1|k+1)=X(k+1|k)+K*(Z(k)-H*X(k+1|k))
    P0=(1-K*H)*P_pred;                  % 误差协方差的更新方程: P(k+1|k+1)=(I-K*P(k+1|k)           
          
end

% 

meassure_err_x = zeros(1,N);
kalman_err_x = zeros(1,N);
kalman_err_v = zeros(1,N);
for k=1:N
    meassure_err_x(k) = Z(k) - X(1,k);
    kalman_err_x(k) = Xkf(1,k) - X(1,k);
    kalman_err_v(k) = sqrt(0.5*(Xkf(2,k) - X(2,k))^2);
end

figure
hold on,box on;
plot(meassure_err_x,'-r.');
plot(kalman_err_x,'-g.');
legend('测量位置','Kalman估计位置');
xlabel('采样时间/s');
ylabel('位置偏差/m');
title('Kalman滤波跟踪')
figure
plot(kalman_err_v);
xlabel('采样时间/s');
ylabel('速度偏差');
title('速度RMSE')
figure
hold on,box on;
plot(Xkf(1,:),'-r.');
plot(X(1,:),'-go');
legend('测量位置','Kalman估计位置');
xlabel('采样时间/s');
ylabel('位置信息/m');
title('Kalman滤波跟踪位置')
figure
hold on,box on;
plot(Xkf(2,:),'-r.');
plot(X(2,:),'-go');
legend('测量速度','Kalman估计位置');
xlabel('采样时间/s');
ylabel('速度信息/m');
title('Kalman滤波跟踪速度')
