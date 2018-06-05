%%体验什么是蒙特卡洛方法 code转自博客http://blog.sciencenet.cn/blog-316653-375888.html
%%蒙特卡洛方法使用随机数的统计特性来模拟数据分布。当所求问题的解是某个事件的概率，或者是某个
%%随机变量的数学期望，或者是与之有关的量时，通过某种试验的方法，得出该事件发生的概率，再通过它得到问题的解。

%%rand(m,n)生成的数据均为(0,1)之间

x=0:0.01:1;y=x.^2;plot(x,y);
Outcome=zeros(4);
staus=10;
for i=1:4  %4次模拟
point=staus.^i; %模拟的随机点数
RandData=rand(2,point); %根据随机点数，产生随机的(x,y)散点,不明白可以试试   
% scatter(RandData(1,:),RandData(2,:))
Below=find(RandData(1,:).^2>RandData(2,:));%寻找位于曲线下的散点
Outcome(i)=length(Below)/length(RandData);%最终结果的表示
end
BelowData=RandData(:,Below);
hold on
scatter(BelowData(1,:),BelowData(2,:))
