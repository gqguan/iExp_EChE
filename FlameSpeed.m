function outTab = FlameSpeed(heightCalibration,firePort,flameHeight)
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
u.La = flameHeight.xm.La*flameHeight.s.La/sqrt(3); % B类不确定度评定，单位：m3/s
u.Lg = flameHeight.xm.Lg*flameHeight.s.Lg/sqrt(3); % B类不确定度评定，单位：m3/s
u.D = sqrt(std(firePort.D)^2+(firePort.a/sqrt(3))^2); % A和B类不确定度评定，单位：m
% 各测量量的自由度
nu.La = 1;
nu.Lg = 1;
nu.D = length(firePort.D);

%% 一元线性回归火焰高度H
n = length(knownX);
nu.H = n-2; % 自由度.火焰高度
lxx = sum((knownX-mean(knownX)).^2);
lxy = sum((knownX-mean(knownX)).*(knownY-mean(knownY)));
lyy = sum((knownY-mean(knownY)).^2);
% 回归系数
b = lxy/lxx;
a = mean(knownY)-b*mean(knownX);
% 回归结果
y = @(x)(a+b*x);
% 残余标准差
s = sqrt(sum((knownY-y(knownX)).^2)/nu.H);
% 回归值的扩展不确定度
U95b = tinv(p,nu.H)*(s*sqrt(1/lxx));
U95a = tinv(p,nu.H)*(s*sqrt(1/n+mean(knownX)^2/lxx));
U95y = @(x)(tinv(p,nu.H)*(s*sqrt(1/n+(x-mean(knownX)).^2/lxx)));
y_uc2 = @(x)(a+b*x+U95y(x));
y_uc1 = @(x)(a+b*x-U95y(x));
eqTxt = sprintf('$\\hat{H} = %.02g+%.4fH$',a,b);
% 火焰高度测量的标准不确定度
u.H = @(x)U95y(x);

%% 火焰传播速度
sn = 4*(Lg+La)/(pi*D*sqrt(D^2+4*H^2));
Fs1 = @(X)eval(subs(diff(sn,Lg)*u.Lg,[Lg La D H],X));
Fs2 = @(X)eval(subs(diff(sn,La)*u.La,[Lg La D H],X));
Fs3 = @(X)eval(subs(diff(sn,D)*u.D,[Lg La D H],X));
Fs4 = @(X)eval(subs(diff(sn,H)*u.H(X(4)),[Lg La D H],X));
uc_s = @(X)sqrt(Fs1(X).^2+Fs2(X).^2+Fs3(X).^2+Fs4(X).^2);
nu_eff = @(X)uc_s(X)^4/(u.Lg^4/nu.Lg+u.La^4/nu.La+u.D^4/nu.D+u.H(X(4))^4/nu.H);
U_s = @(X)round(tinv(p,nu_eff(X))*uc_s(X),2,'significant');

%% 燃气浓度
xg = Lg/(Lg+La);
Fx1 = @(X)eval(subs(diff(xg,Lg)*u.Lg,[Lg La D H],X));
Fx2 = @(X)eval(subs(diff(xg,La)*u.La,[Lg La D H],X));
uc_x = @(X)sqrt(Fx1(X).^2+Fx2(X).^2);
U_x = @(X)uc_x(X)*2;

%% 输出1：火焰高度校正工作曲线
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

%% 输出2：不同燃气浓度下的法向火焰传播速度
Sval = zeros(1,length(expH));
Suc = zeros(1,length(expH));
Sn = string;
Xg = string;
xval = zeros(1,length(expH));
F = zeros(4,length(expH));
for i = 1:length(expH)
    X = [expLg(i),expLa(i),expD,expH(i)];
    F(:,i) = [Fs1(X),Fs2(X),Fs3(X),Fs4(X)]';
    Sval(i) = eval(subs(sn,[Lg La D H],X));
    Suc(i) = U_s(X);
    [~,Sn(i)] = FixValue(Sval(i),Suc(i));
    Sn(i) = strcat(Sn(i)," m/s (p=0.95)");
    xval(i) = eval(subs(xg,[Lg La D H],X));
    [~,Xg(i)] = FixValue(xval(i),U_x(X));
%     xval(i) = expLg(i)/(expLg(i)+expLa(i));
end
figure(2)
errorbar(xval,Sval,Suc,'o')
xlabel('$x_\mathrm{g}$','interpreter','latex')
ylabel('$s_\mathrm{n}$ (m/s)','interpreter','latex')
% 列表输出
outTab = table(Xg',Sn','VariableNames',{'xg','sn'});

%% 输出3：各测量量对结果的不确定度贡献
figure(3)
bar(F)
xticklabels({'Lg' 'La' 'D' 'H'})
xlabel('Measured Variables')
ylabel('uc_i')
