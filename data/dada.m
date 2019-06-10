%4Q
N=64;
node=(sqrt(N)+1)^2;
Nodenum=linspace(1,node,node);
XYZ=zeros(node,7);
XYZ(:,1)=Nodenum';
bcn=zeros(1,sqrt(N)+1);
for i=1:(sqrt(N)+1)
    bcn(i)=(sqrt(N)+1)*(i-1)+1;
    XYZ(bcn(i),2:3)=1;
end
XYZ(:,4)=1;
for i=1:(sqrt(N)+1)
    for j=1:(sqrt(N)+1)
        nd=(i-1)*(sqrt(N)+1)+j;
        x=(j-1)*(2/(sqrt(N)));
        y=(i-1)*(2/(sqrt(N)));
        XYZ(nd,5)=x;
        XYZ(nd,6)=y;
    end
end
dlmwrite('XYZ.txt', XYZ, 'delimiter', '\t','precision', 3,'newline', 'pc');

Force=zeros((sqrt(N)+1),3);
for i=1:(sqrt(N)+1)
    Force(i,1)=i*(sqrt(N)+1);
    Force(i,2)=1;
end
Force(:,3)=20*(2/(sqrt(N)));
Force(1,3)=20*(2/(sqrt(N)))/2;
Force((sqrt(N)+1),3)=20*(2/(sqrt(N)))/2;
dlmwrite('Force.txt', Force, 'delimiter', '\t','precision', 3,'newline', 'pc');

LM=zeros(N,6);
LM(:,6)=1;
LM(:,1)=linspace(1,N,N);
for i=1:sqrt(N)
    for j=1:sqrt(N)
        node1=(i-1)*(sqrt(N)+1)+j;
        node2=(i-1)*(sqrt(N)+1)+j+1;
        node3=i*(sqrt(N)+1)+j+1;
        node4=i*(sqrt(N)+1)+j;
        ne=(i-1)*sqrt(N)+j;
        LM(ne,2)=node1;
        LM(ne,3)=node2;
        LM(ne,4)=node3;
        LM(ne,5)=node4;
    end
end
dlmwrite('LM.txt', LM, 'delimiter', '\t','precision', 9,'newline', 'pc');

%Plate
N=64;
node=(sqrt(N)+1)^2;
Nodenum=linspace(1,node,node);
XYZ=zeros(node,10);
XYZ(:,1)=Nodenum';
bcn=zeros(1,sqrt(N)+1);
XYZ(1,4:6)=1;
XYZ(sqrt(N)+1,4:6)=1;
XYZ(((sqrt(N)+1)^2-sqrt(N)),4:6)=1;
XYZ((sqrt(N)+1)^2,4:6)=1;

XYZ(:,2)=1;
XYZ(:,3)=1;
XYZ(:,7)=1;
for i=1:(sqrt(N)+1)
    for j=1:(sqrt(N)+1)
        nd=(i-1)*(sqrt(N)+1)+j;
        x=(j-1)*(2/(sqrt(N)));
        y=(i-1)*(2/(sqrt(N)));
        XYZ(nd,8)=x;
        XYZ(nd,9)=y;
    end
end
dlmwrite('XYZ.txt', XYZ, 'delimiter', '\t','precision', 9,'newline', 'pc');


LM=zeros(N,6);
LM(:,6)=1;
LM(:,1)=linspace(1,N,N);
for i=1:sqrt(N)
    for j=1:sqrt(N)
        node1=(i-1)*(sqrt(N)+1)+j;
        node2=(i-1)*(sqrt(N)+1)+j+1;
        node3=i*(sqrt(N)+1)+j+1;
        node4=i*(sqrt(N)+1)+j;
        ne=(i-1)*sqrt(N)+j;
        LM(ne,2)=node1;
        LM(ne,3)=node2;
        LM(ne,4)=node3;
        LM(ne,5)=node4;
    end
end
dlmwrite('LM.txt', LM, 'delimiter', '\t','precision', 3,'newline', 'pc');




syms x0
syms y0
syms a
syms b
syms mu
syms x
syms y
D=[1 mu 0;mu 1 0;0 0 (1-mu)/2];
Q=[0 0 0 2 0 0 6*x 2*y 0 0 6*x*y 0;0 0 0 0 0 2 0 0 2*x 6*y 0 6*x*y;0 0 0 0 2 0 0 4*x 4*y 0 6*x.^2 6*y.^2];
J=Q.'*D*Q;
J1=int(int(J,x,x0-a,x0+a),y,y0-b,y0+b);
M(1,:)=[1,(x0-a/2),(y0-b/2),(x0-a/2).^2,(x0-a/2).*(y0-b/2),(y0-b/2).^2,(x0-a/2).^3,(x0-a/2).^2.*(y0-b/2),(x0-a/2).*(y0-b/2).^2,(y0-b/2).^3,(x0-a/2).^3.*(y0-b/2),(x0-a/2).*(y0-b/2).^3];
M(2,:)=[0,0,1,0,(x0-a/2),2*(y0-b/2),0,(x0-a/2).^2,2*(x0-a/2)*(y0-b/2),3*(y0-b/2).^2,(x0-a/2).^3,3*(x0-a/2).*(y0-b/2).^2];
M(3,:)=[0,-1,0,-2*(x0-a/2),-(y0-b/2),0,-3*(x0-a/2).^2,-2*(x0-a/2)*(y0-b/2),-(y0-b/2).^2,0,-3*(x0-a/2).^2.*(y0-b/2),-(y0-b/2).^3];
M(4,:)=[1,(x0+a/2),(y0-b/2),(x0-a/2).^2,(x0+a/2).*(y0-b/2),(y0-b/2).^2,(x0+a/2).^3,(x0+a/2).^2.*(y0-b/2),(x0+a/2).*(y0-b/2).^2,(y0-b/2).^3,(x0+a/2).^3.*(y0-b/2),(x0+a/2).*(y0-b/2).^3];
M(5,:)=[0,0,1,0,(x0+a/2),2*(y0-b/2),0,(x0+a/2).^2,2*(x0+a/2)*(y0-b/2),3*(y0-b/2).^2,(x0+a/2).^3,3*(x0+a/2).*(y0-b/2).^2];
M(6,:)=[0,-1,0,-2*(x0+a/2),-(y0-b/2),0,-3*(x0+a/2).^2,-2*(x0+a/2)*(y0-b/2),-(y0-b/2).^2,0,-3*(x0+a/2).^2.*(y0-b/2),-(y0-b/2).^3];
M(7,:)=[1,(x0+a/2),(y0+b/2),(x0+a/2).^2,(x0+a/2).*(y0+b/2),(y0+b/2).^2,(x0+a/2).^3,(x0+a/2).^2.*(y0+b/2),(x0+a/2).*(y0+b/2).^2,(y0+b/2).^3,(x0+a/2).^3.*(y0+b/2),(x0+a/2).*(y0+b/2).^3];
M(8,:)=[0,0,1,0,(x0+a/2),2*(y0+b/2),0,(x0+a/2).^2,2*(x0+a/2)*(y0+b/2),3*(y0+b/2).^2,(x0+a/2).^3,3*(x0+a/2).*(y0+b/2).^2];
M(9,:)=[0,-1,0,-2*(x0+a/2),-(y0+b/2),0,-3*(x0+a/2).^2,-2*(x0+a/2)*(y0+b/2),-(y0+b/2).^2,0,-3*(x0+a/2).^2.*(y0+b/2),-(y0+b/2).^3];
M(10,:)=[1,(x0-a/2),(y0+b/2),(x0-a/2).^2,(x0-a/2).*(y0+b/2),(y0+b/2).^2,(x0-a/2).^3,(x0-a/2).^2.*(y0+b/2),(x0-a/2).*(y0+b/2).^2,(y0+b/2).^3,(x0-a/2).^3.*(y0+b/2),(x0-a/2).*(y0+b/2).^3];
M(11,:)=[0,0,1,0,(x0-a/2),2*(y0+b/2),0,(x0-a/2).^2,2*(x0-a/2)*(y0+b/2),3*(y0+b/2).^2,(x0-a/2).^3,3*(x0-a/2).*(y0+b/2).^2];
M(12,:)=[0,-1,0,-2*(x0-a/2),-(y0+b/2),0,-3*(x0-a/2).^2,-2*(x0-a/2)*(y0+b/2),-(y0+b/2).^2,0,-3*(x0-a/2).^2.*(y0+b/2),-(y0+b/2).^3];
K=(inv(M)).'*J1*inv(M);



K1=zeros(78,4);
K1(:,1)=linspace(0,77,78);
for i=1:12
    for j=1:i
        nd=i*(i-1)/2+(i-j)+1;
        K1(nd,2)=j;
        K1(nd,3)=i;
    end
end
LM=[3 4 5 9 10 11 15 16 17 21 22 23];
for l=1:78
    K1(l,2)=LM(K1(l,2));
    K1(l,3)=LM(K1(l,3));
    K1(l,4)=K1(l,3)*(K1(l,3)-1)/2+(K1(l,3)-K1(l,2));
end


%Subpara
N=16;
node=(2*sqrt(N)+1)^2;
Nodenum=linspace(1,node,node);
XYZ=zeros(node,7);
XYZ(:,1)=Nodenum';
bcn=zeros(1,2*sqrt(N)+1);
for i=1:(sqrt(N)+1)
    bcn(i)=(sqrt(N)+1)*(i-1)+1;
    XYZ(bcn(i),2:3)=1;
end
XYZ(:,4)=1;
for i=1:(sqrt(N)+1)
    for j=1:(sqrt(N)+1)
        nd=(i-1)*(sqrt(N)+1)+j;
        x=(j-1)*(2/(sqrt(N)));
        y=(i-1)*(2/(sqrt(N)));
        XYZ(nd,5)=x;
        XYZ(nd,6)=y;
    end
end
dlmwrite('XYZ.txt', XYZ, 'delimiter', '\t','precision', 3,'newline', 'pc');

Force=zeros((sqrt(N)+1),3);
for i=1:(sqrt(N)+1)
    Force(i,1)=i*(sqrt(N)+1);
    Force(i,2)=1;
end
Force(:,3)=20*(2/(sqrt(N)));
Force(1,3)=20*(2/(sqrt(N)))/2;
Force((sqrt(N)+1),3)=20*(2/(sqrt(N)))/2;
dlmwrite('Force.txt', Force, 'delimiter', '\t','precision', 3,'newline', 'pc');

LM=zeros(N,6);
LM(:,6)=1;
LM(:,1)=linspace(1,N,N);
for i=1:sqrt(N)
    for j=1:sqrt(N)
        node1=(i-1)*(sqrt(N)+1)+j;
        node2=(i-1)*(sqrt(N)+1)+j+1;
        node3=i*(sqrt(N)+1)+j+1;
        node4=i*(sqrt(N)+1)+j;
        ne=(i-1)*sqrt(N)+j;
        LM(ne,2)=node1;
        LM(ne,3)=node2;
        LM(ne,4)=node3;
        LM(ne,5)=node4;
    end
end
dlmwrite('LM.txt', LM, 'delimiter', '\t','precision', 3,'newline', 'pc');