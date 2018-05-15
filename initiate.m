function [ T ] = initiate(n,m,t0 )
T=zeros(m,n);
for i=1:n
    for j=1:m
        T(i,j)=t0;
    end   
end
end

