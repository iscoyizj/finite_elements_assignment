% ’¡≤¬ 
clear;
XYZ=textread('XYZ.txt');
LM=textread('LM.txt');
n1=size(XYZ);
n2=size(LM);
F=zeros(n1(1),3);
F(:,1)=linspace(1,n1(1),n1(1));
F(:,2)=1;
a=9;
for i=1:n2(1)
    x1=XYZ(LM(i,2),5);
    x2=XYZ(LM(i,3),5);
    h=x2-x1;
    F(LM(i,2),3)=F(LM(i,2),3)+h^2/4*((x1+x2)/2*a/9-a*h/18);
    F(LM(i,3),3)= F(LM(i,3),3)+h^2/4*((x1+x2)/2*a/9+a*h/18);
    F(LM(i,4),3)= F(LM(i,4),3)+h^2/4*((x1+x2)/2*a/9+a*h/18);
    F(LM(i,5),3)= F(LM(i,5),3)+h^2/4*((x1+x2)/2*a/9-a*h/18);
    F(LM(i,6),3)= F(LM(i,6),3)+h^2/4*((x1+x2)/2*4*a/9);
    F(LM(i,7),3)= F(LM(i,7),3)+h^2/4*((x1+x2)/2*4*a/9+2*a*h/9);
    F(LM(i,8),3)= F(LM(i,8),3)+h^2/4*((x1+x2)/2*4*a/9);
    F(LM(i,9),3)= F(LM(i,9),3)+h^2/4*((x1+x2)/2*4*a/9-2*a*h/9);
    F(LM(i,10),3)= F(LM(i,10),3)+h^2/4*((x1+x2)/2*16*a/9);
end
dlmwrite('Force.txt', F, 'delimiter', '\t','precision', 3,'newline', 'pc');

%∑÷∆¨≤‚ ‘
clear;
XYZ=textread('XYZ.txt');
LM=textread('LM.txt');
n1=size(XYZ);
n2=size(LM);
F=zeros(n1(1),3);
F(:,1)=linspace(1,n1(1),n1(1));
F(:,2)=1;
a=9;
for i=1:n2(1)
    x1=XYZ(LM(i,2),5);
    x2=XYZ(LM(i,3),5);
    h=x2-x1;
    F(LM(i,2),3)=F(LM(i,2),3)+a*h^2/36;
    F(LM(i,3),3)= F(LM(i,3),3)+a*h^2/36;
    F(LM(i,4),3)= F(LM(i,4),3)+a*h^2/36;
    F(LM(i,5),3)= F(LM(i,5),3)+a*h^2/36;
    F(LM(i,6),3)= F(LM(i,6),3)+a*h^2/9;
    F(LM(i,7),3)= F(LM(i,7),3)+a*h^2/9;
    F(LM(i,8),3)= F(LM(i,8),3)+a*h^2/9;
    F(LM(i,9),3)= F(LM(i,9),3)+a*h^2/9;
    F(LM(i,10),3)= F(LM(i,10),3)+4*a*h^2/9;
end
dlmwrite('Force.txt', F, 'delimiter', '\t','precision', 3,'newline', 'pc');
