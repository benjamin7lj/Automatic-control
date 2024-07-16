clc;clear ; close all; 
syms x1 x2 x3 x4 u
m=0.125; L=16.75; r=21.5; g=9.81; Rm=2.6; Km=0.0076; Kt=0.0076; Kg=70;
Beq=0.004; nm=0.87; ng=0.85; Vm=6; Jeq=0.0023;
%Tout=nm*ng*Kt*Kg*(Vm-Kg*Km*x1)/Rm;
Tout=u;
f1=x2; f3=x4;
f2=(Tout-Beq*x2+(3/4)*m*r*g*x3)/(Jeq+m*r^2/4);
f4=(Tout-Beq*x2+(Jeq+m*r^2)*(g/r)*x3)/((Jeq+m*r^2)*(4/3)*(L/r)-m*L*r);
Jx=jacobian([f1,f2,f3,f4],[x1,x2,x3,x4]);
Ju=jacobian([f1,f2,f3,f4],u);
A=subs(Jx,[x1,x2,x3,x4],[0,0,0,0]);
B=subs(Ju,u,0);
[V,D]=eig(A);
lambda=eig(A);