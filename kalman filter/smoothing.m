% Kalman
clear;
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
S=zeros(N,2);
SN=[0,0;0,0];
Xkf=zeros(2,N);
Xkf(:,1)=X(:,1);

Q = [0,1;1,0];
R =1;
W = sqrt(Q)*randn(2,N);                               % 模拟产生过程噪声(高斯分布的随机噪声)
V = sqrt(R)*randn(1,N);                               % 模拟产生测量噪声
S=zeros(N,2);

for k=2:N

    X(:,k) = A * X(:,k-1)+g*B+V(k-1);               % 状态方程:X(k+1)=Φ(k)X(k)+G(k)V(k),其中G(k)=1
    Z(:,k)=H*X(:,k)+V(k);                           % 观测方程:Z(k+1)=H(k+1)X(k+1)+W(k+1),Z(k+1)是k+1时刻的观测值
end
err_p(1,1)=P0(1,1);  
err_p(1,2)=P0(2,2);   
%% 前向滤波器
for k=2:N
    
    X_pred=A*Xkf(:,k-1)+g*B;
    P_pred=A*P0*A'+Q;                   % 一步预测的协方差 P(k+1|k)
    K=P_pred*H'*inv(H*P_pred*H'+R);     % 卡尔曼滤波器增益 K
    Xkf(:,k)=X_pred+K*(Z(k)-H*X_pred);  % 状态更新方程 X(k+1|k+1)=X(k+1|k)+K*(Z(k)-H*X(k+1|k))
    P0=(1-K*H)*P_pred;                  % 误差协方差的更新方程: P(k+1|k+1)=(I-K*P(k+1|k)           
    err_p(k,1)=P0(1,1);                 %位移误差均方值  
    err_p(k,2)=P0(2,2);                 %速度误差均方值
      
end
z_hatN=0;
z_hat=zeros(1,N);
z_hat(N,:)=z_hatN;
%% 反向
for k=1:N
    S_back=S(:,N+1-k)+H'*R^-1*H;
    z=z_hat(N+1-k)+H'*R^-1*Z(:,N+1-k);
    K=S_back*(S_back+Q^-1)^-1;
    S=A'*(1-K)*S*(1-K)'*A+A'*K*Q^-1*K*A;
    z_hat=A'*(1-K)*z;
end
PP=[0,0;0,0];
%% 平滑计算；
for k=1:N
    PP=(P0^-1+S)^-1;
    X(:,N-k)=P*(P0^-1*Xkf+z_hat);
end

plot(t,P(1:N),'k-+',t,PP(1:N),'r-x',t,P0,'b-*');
xlabel('步长');
ylabel('状态');
legend('估计方差','平滑方差','一步预测方差');
grid on
