%% 例题9-1
% 试列出实验数据的线性拟合预测结果
% by Dr. Guan Guoqiang @SCUT on 2021/9/3

%% 初始化
clear
knownX = [0.602 4.250 3.550 1.432]'/100; % 单位：m
knownY = [0.800 4.400 3.600 1.600]'/100; % 单位：m
p = 0.95;
syms Lg La D H
expH = [0.01 0.0085 0.01218 0.0142 0.0198 0.02042 0.0253 0.028]; % 单位：m
expLa = [0.175 0.14 0.17 0.235 0.26 0.315 0.34 0.375]/3600; % 单位：m3/s
expLg = [0.0045 0.0042 0.0054 0.0066 0.0078 0.009 0.0102 0.0114]/3600; % 单位：m3/s
expD = 4.15e-3; % 单位：m
u.La = 0.6*0.04/3600/sqrt(3); % 单位：m3/s
u.Lg = 250e-6/60/sqrt(3); % 单位：m3/s
u.D = 0.03e-3/sqrt(3); % 单位：m
% 各测量量的自由度
nu.La = 1;
nu.Lg = 1;
nu.D = 4;

%% 一元线性回归火焰高度H
n = length(knownX);
nu.H = n-2; % 自由度.火焰高度
lxx = sum((knownX-mean(knownX)).^2);
lxy = sum((knownX-mean(knownX)).*(knownY-mean(knownY)));
lyy = sum((knownY-mean(knownY)).^2);
% 回归系数
b = lxy/lxx;
a = mean(knownY)-b*mean(knownX);
% 残余标准差
s = sqrt((lyy-lxy^2/lxx)/nu.H);
% 回归值的扩展不确定度
U95b = tinv(p,nu.H)*(s*sqrt(1/lxx));
U95a = tinv(p,nu.H)*(s*sqrt(1/n+mean(knownX)^2/lxx));
U95y = @(x)(tinv(p,nu.H)*(s*sqrt(1/n+(x-mean(knownX)).^2/lxx)));
% 回归结果
y = @(x)(a+b*x);
y_uc2 = @(x)(a+b*x+U95y(x));
y_uc1 = @(x)(a+b*x-U95y(x));
eqTxt = sprintf('$\\hat{H} = %.02g+%.4fH$',a,b);
% 火焰高度测量的标准不确定度
u.H = @(x)(y_uc2(x)-y_uc1(x));

%% 火焰传播速度
s = 4*(Lg+La)/(pi*D*sqrt(D^2+4*H^2));
F1 = @(X)eval(subs(diff(s,Lg)*u.Lg,[Lg La D H],X));
F2 = @(X)eval(subs(diff(s,La)*u.La,[Lg La D H],X));
F3 = @(X)eval(subs(diff(s,D)*u.D,[Lg La D H],X));
F4 = @(X)eval(subs(diff(s,H)*u.H(X(4)),[Lg La D H],X));
uc = @(X)sqrt(F1(X).^2+F2(X).^2+F3(X).^2+F4(X).^2);
nu_eff = @(X)uc(X)^4/(u.Lg^4/nu.Lg+u.La^4/nu.La+u.D^4/nu.D+u.H(X(4))^4/nu.H);
U = @(X)tinv(p,nu_eff(X))*uc(X);

%% 输出1
xmesh = linspace(min(knownX),max(knownX));
figure(1)
hold on
plot(xmesh,y(xmesh),'b')
plot(knownX,knownY,'bo')
plot(xmesh,y_uc1(xmesh),'r--')
plot(xmesh,y_uc2(xmesh),'r--')
xlabel('$H$ (m)','interpreter','latex')
ylabel('$\hat{H}$ (m)','interpreter','latex')
text(min(xmesh)*4,min(y(xmesh))*3,eqTxt,'interpreter','latex')
hold off

%% 输出2
Sval = zeros(1,length(expH));
Suc = zeros(1,length(expH));
xval = zeros(1,length(expH));
for i = 1:length(expH)
    X = [expLg(i),expLa(i),expD,expH(i)];
    Sval(i) = eval(subs(s,[Lg La D H],X));
    Suc(i) = U(X);
    xval(i) = expLg(i)/(expLg(i)+expLa(i));
end
figure(2)
errorbar(xval,Sval,Suc,'o')
xlabel('$x_\mathrm{g}$','interpreter','latex')
ylabel('$s_\mathrm{n}$','interpreter','latex')