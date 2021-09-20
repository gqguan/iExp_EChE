%% 例题9-1
% 试列出实验数据的线性拟合预测结果
% by Dr. Guan Guoqiang @SCUT on 2021/9/3

%% 初始化
knownX = [180 123 141 108 121 180 145 191 190 155 165 143 116 134 110 ...
    150 107 147 115 144 153 161 177 104 151 204 158 151 127 141 190 190 ...
    177 154]';
knownY = [200 110 125 110 125 240 165 205 190 160 195 160 100 135 130 ...
    170 115 155 120 160 145 145 205 100 180 235 130 135 135 135 220 210 ...
    185 150]';
p = 0.95;

%% 一元线性回归
n = length(knownX);
nu = n-2; % 自由度
lxx = sum((knownX-mean(knownX)).^2);
lxy = sum((knownX-mean(knownX)).*(knownY-mean(knownY)));
lyy = sum((knownY-mean(knownY)).^2);
% 回归系数
b = lxy/lxx;
a = mean(knownY)-b*mean(knownX);
% 残余标准差
s = sqrt((lyy-lxy^2/lxx)/nu);
% 回归值的扩展不确定度
U95b = tinv(p,nu)*(s*sqrt(1/lxx));
U95a = tinv(p,nu)*(s*sqrt(1/n+mean(knownX)^2/lxx));
U95y = @(x)(tinv(p,nu)*(s*sqrt(1/n+(x-mean(knownX)).^2/lxx)));
% 回归结果
y = @(x)(a+b*x);
y_uc2 = @(x)(a+b*x+U95y(x));
y_uc1 = @(x)(a+b*x-U95y(x));
eqTxt = sprintf('y = %.2f+%.2fx',a,b);

%% 输出
xmesh = linspace(min(knownX),max(knownX));
figure(1)
hold on
plot(xmesh,y(xmesh),'b')
plot(xmesh,y_uc1(xmesh),'r--')
plot(xmesh,y_uc2(xmesh),'r--')
xlabel('$x_i$','interpreter','latex')
ylabel('$y_i$','interpreter','latex')
text(min(xmesh)*1.8,min(y(xmesh))*1.8,eqTxt)