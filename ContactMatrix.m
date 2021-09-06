load('ContactMatrix16');

cinter=zeros(16,4);
for i=1:16
    cinter(i,1)=sum(C16(i,1:4));
    cinter(i,2)=sum(C16(i,5:9));
    cinter(i,3)=sum(C16(i,10:13));
    cinter(i,4)=sum(C16(i,14:16));
end
c4=zeros(4,4);
for i=1:4
    c4(1,i)=sum(cinter(1:4,i));
    c4(2,i)=sum(cinter(5:9,i));
    c4(3,i)=sum(cinter(10:13,i));
    c4(4,i)=sum(cinter(14:16,i));
end