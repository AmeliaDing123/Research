function [ dist] = Euclidean(A,B )
%E �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
N=size(B,1);
for i=1:N
    S(i,:)=(A-B(i,:)).^2;
end

dist=sqrt(sum(S,2));


end

