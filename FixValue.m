%% 按输入的不确定度圆整输入数值
%
% by Dr. Guan Guoqiang @ SCUT on 2021-9-22
function [valOut,strOut] = FixValue(valIn,valU)
% 将不确定度结果转为字符
strUc = sprintf('%.2g',valU);
% 字符串strUc中最后第2个字符若为0则在strUc中再补一个有效数字0
if isequal(strUc(end-1),'0')
    strUc = [strUc,'0'];
end
strUcFra = strUc(strfind(strUc,'.')+1:end);
% 不确定度的有效数字
sigNum = length(strUc)-strfind(strUc,'.');
% 将输入数值转换为字符
fmt1 = sprintf('%%.%df',sigNum+1);
strVal0 = num2str(valIn,fmt1);
strVal0Int = strVal0(1:strfind(strVal0,'.')-1);
strVal0Fra = strVal0(strfind(strVal0,'.')+1:end);
% valueIn的小数长度应不小于valueUc的长度+1
if length(strVal0Fra) >= length(strUcFra)+1
    strVal = [strVal0Int,'.',strVal0Fra(1:length(strUcFra)+1)];
else
    warning('输入数值长度小于不确定度声明！')
    return
end
% 按“四舍六入五留双”规则取整
fmt2 = sprintf('%%0%dd',sigNum);
valueOut = sprintf(fmt2,convergent(str2double(strVal)*10^sigNum));
if isempty(valueOut(1:end-sigNum))
    strOut = ['0.',valueOut(end-sigNum+1:end)];
else
    strOut = [valueOut(1:end-sigNum),'.',valueOut(end-sigNum+1:end)];
end
valOut = str2double(strOut);
strOut = [strOut,'±',strUc];
end