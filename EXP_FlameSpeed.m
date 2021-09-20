%% 燃气法向火焰传播速度测定实验
%
% by Dr. Guan Guoqiang @ SCUT on 2021-9-18

%% 初始化
clear
flameHeight.s.Lg = 0.04; % 燃气转子流量计精度等级
flameHeight.s.La = 0.04; % 空气转子流量计精度等级
flameHeight.xm.Lg = SIConvert(1.6,'LPM'); % 燃气转子流量计量程
flameHeight.xm.La = SIConvert(1.0,'m3/h'); % 空气转子流量计量程
firePort.a = SIConvert(0.03,'mm'); % 游标卡尺最大允许误差

%% 实验结果
% 测高仪测量校正实验
heightCalibration.H = SIConvert([4.412 0.804 1.604 2.502 4.810]','cm');
heightCalibration.H1 = SIConvert([25.850 22.380 23.010 24.010 25.090]','cm');
heightCalibration.H2 = SIConvert([21.492 21.498 21.498 21.524 20.322]','cm');
% 火孔内径
firePort.D = SIConvert([1.350 1.344 1.334 1.344],'cm');
% 火焰高度测定
flameHeight.Lg = SIConvert([0.28 0.36 0.44 0.48 0.56 0.72 0.80]','LPM');
flameHeight.La = SIConvert([0.50 0.50 0.55 0.65 0.75 0.85 1.00]','m3/h');
flameHeight.H1 = SIConvert([23.490 23.246 23.594 23.930 24.350 24.848 25.344]','cm');
flameHeight.H2 = SIConvert([21.550 21.538 21.512 21.508 21.510 21.510 21.510]','cm');

%% 计算火焰法向传播速度并绘图输出结果
FlameSpeed(heightCalibration,firePort,flameHeight)