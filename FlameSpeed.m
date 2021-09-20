function FlameSpeed(heightCalibration,firePort,flameHeight)
%% 关联测高仪工作曲线、计算法向火焰传播速度和不确定度
% （实验记录表Ref.#210905R）
% 
% by Dr. Guan Guoqiang @SCUT on 2021/9/3

%% 初始化
knownX = abs(heightCalibration.H1-heightCalibration.H2); % 单位：m
knownY = heightCalibration.H; % 单位：m
p = 0.95;
syms Lg La D H
expH = abs(flameHeight.H1-flameHeight.H2); % 单位：m
expLa = flameHeight.La; % 单位：m3/s
expLg = flameHeight.Lg; % 单位：m3/s
expD = mean(firePort.D); % 单位：m
u.La = flameHeight.xm.La*flameHeight.s.La/sqrt(3); % 单位：m3/s
u.Lg = flameHeight.xm.Lg*flameHeight.s.Lg/sqrt(3); % 单位：m3/s
u.D = firePort.a/sqrt(3); % 单位：m
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

%% 输出3：各测量量对结果的不确定度贡献
figure(3)
bar([F1(X),F2(X),F3(X),F4(X)])
xticklabels({'Lg' 'La' 'D' 'H'})
