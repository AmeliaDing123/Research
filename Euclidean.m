function [ dist] = Euclidean(A,B )
%E 此处显示有关此函数的摘要
%   此处显示详细说明
N=size(B,1);
for i=1:N
    S(i,:)=(A-B(i,:)).^2;
end

dist=sqrt(sum(S,2));


end

