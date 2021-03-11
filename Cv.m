function [ COVMAT ] = Cv( inputData,E )
% ? covarianceMatrix( inputData )
% ? 这是一个计算协方差矩阵的函数
% ? inputData ? 输入数据
% ? 每一行为一个维度
% ? 每一列为一个样本
%获得输入数据维度
[m,n] = size(inputData');
%创建协方差矩阵
COVMAT = zeros(m,m);

%取得每维数据平均值

%计算协方差
for i = 1:n
    
        
            COVMAT= ((inputData(i,:)-E)'*(inputData(i,:)-E))+COVMAT;
        
    
end
COVMAT=COVMAT/n;
