function [ COVMAT ] = Cv( inputData,E )
% ? covarianceMatrix( inputData )
% ? ����һ������Э�������ĺ���
% ? inputData ? ��������
% ? ÿһ��Ϊһ��ά��
% ? ÿһ��Ϊһ������
%�����������ά��
[m,n] = size(inputData');
%����Э�������
COVMAT = zeros(m,m);

%ȡ��ÿά����ƽ��ֵ

%����Э����
for i = 1:n
    
        
            COVMAT= ((inputData(i,:)-E)'*(inputData(i,:)-E))+COVMAT;
        
    
end
COVMAT=COVMAT/n;
