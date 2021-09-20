%% 例题6-4
% 用弓高弦长法测量大工件。车间工人用卡尺量得h=50mm、l=500mm；质检部门用高精度
% 等级的卡尺测量的h=50.1mm、l=499mm。问车间工人测量该工件的系统误差及修正后的
% 测量结果。
%
% by Dr. Guan Guoqiang @SCUT on 2021/9/2

%% 初始化
clear
syms x y
l1 = 500; h1 = 50;
lc = 499; hc = 50.1;
dh = hc-h1;
dl = lc-l1;

%% 
% 定义直径计算式
D = @(l,h)(l^2/4/h+h);
% 定义直径偏差的计算式
dD = @(l,h)(eval(subs(diff(D(x,y),x)*dl+diff(D(x,y),y)*dh,[x,y],[l h])));
% 修正直径结果
cD = @(l,h)(D(l,h)+dD(l,h));

%% 输出结果
fprintf('直径测量的系统误差为%.1fmm，修正后测量直径为%.1fmm\n',dD(l1,h1),...
    cD(l1,h1))

%% 例题6-5
% 已知弓高的标准差为0.005mm、弦长的标准差为0.01mm，求工件直径的标准差和修正后
% 的测量结果
sigma.h = 0.005; sigma.l = 0.01;
sigma.D = @(l,h)(eval(sqrt(subs(diff(D(x,y),x)^2*sigma.l^2 ...
    +diff(D(x,y),y)^2*sigma.h^2,[x,y],[l h]))));

%% 输出结果
fprintf('修正后直径测量的标准差为%.2fmm\n',sigma.D(l1,h1))


