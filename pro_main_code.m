clc;
clear;
close all;
A=[0 1 0 0;0 -15.47 41.68 0;0 0 0 1;0 -14.89 84.05 0];
B=[0;27.12;0;23.13];
C=[1 0 0 0];
D=0;

Q=ctrb(A,B);
N=obsv(A,C);
% Eigen value % vector
[V,E]=eig(A);

%jordan form
V1=inv(V);
An=V1*A*V;
Bn=V1*B;
Cn=C*V;
%cnvert state space to transfer fc
ltiSys=ss(A,B,C,D);
SYS=tf(ltiSys);

%poles and zeros plot
h=pzplot(SYS);
grid on;
%bode diagram
figure()
bode(SYS)
%nyquist diagram
figure()
nyquist(SYS)
%lyapunov
% Define the Lyapunov function
Ql = eye(4); % Positive-definite matrix

% Solve the Lyapunov equation
X = lyap(A', Ql);

% Check stability
if all(eig(X) > 0)
    disp('The system is stable.');
else
    disp('The system is not stable.');
end

%transmision state matrix
syms s t
PHI=inv(s*eye(4)-A);
phi=ilaplace(PHI);

phi_t=expm(t*A);

X_0=[0;0;1;0];
u=[1/s];
state_vec=ilaplace(PHI*X_0)+ilaplace(PHI*B*u);
y=C*state_vec;

%feedback

des_pole=[-5;-3;-20;-10];
K_acker=acker(A,B,des_pole);

%S_info=stepinfo(out.data,out.time);
%static compensator
Gcl0=-C*inv((A-B*K_acker))*B;
Gcl0_inv=inv(Gcl0);

%integral feedback
A1=[A zeros(4,1);-C 0];
B1=[B;0];
C1=[C 0];
des_pole1=[-5;-3;-10;-8;-6];
K_acker1=acker(A1,B1,des_pole1);
K_11=K_acker1(1,1:4);
K_12=K_acker1(1,5);
%zero in 0
A2=[0 1 0 0;0 0 1 0;0 0 0 1;0 679.6 84.05 15.47];
B2=[0;0;0;1];
C2=[0 -1315 0 27.12];
A12=[A2 zeros(4,1);-C 0];
B12=[B2;0];
C12=[C2 0];
des_pole1=[-50;-30;-100;-80;-60];
K_acker12=acker(A12,B12,des_pole1);

%param change
x=0.2;
Ap=A+x*A;
Bp=B+x*B;
Cp=C+x*C;



%observer
des_pole_obs=[-20;-3;-80;-40];
l=acker(A',C',des_pole_obs);
L=l';

des_pole_obs1=[-15;-8;-25;-20];
l1=acker(A',C',des_pole_obs1);
L1=l1';

%lqr-feedback

QQ=[10 0 0 0;0 1 0 0;0 0 100 0;0 0 0 1]; R=10;
K_lqr=lqr(A,B,QQ,R);

%static compensator
Gcl0_lqr=-C*inv((A-B*K_lqr))*B;
Gcl0_inv_lqr=inv(Gcl0);

%integral feedback
QQ1=[10 0 0 0 0;0 1 0 0 0;0 0 100 0 0;0 0 0 1 0;0 0 0 0 10]; 

K_lqr1=lqr(A1,B1,QQ1,R);
K_11_lqr=K_lqr1(1,1:4);
K_12_lqr=K_lqr1(1,5);


%lqr_observer
Ql=[10 0 0 0;0 1 0 0;0 0 100 0;0 0 0 1]; Rl=1;
l_lqr=lqr(A',C',Ql,Rl);
L_lqr=l_lqr';

%minimum_order_observer
%DD=[-10 0;0 -10];
%Q_MOO=[0 0 1 0;1 0 0 0;0 0 0 1;0 1 0 0];
%A_MOO=inv(Q)*A*Q;
%A11_M=[0 0;1 0];A21_M=[0 1;0 0];A12_M=[0 0;0 679.6383];A22_M=[0 8405;1 -1547];

%L_Moo_hat=(DD-A11_M)*inv(A21_M);
%L_MOO=[1 0 0 0;0 1 0 0;0 0 1 -10]*inv(Q);
%R_MOO=L_MOO*B;
%F_MOO=inv([C;L_MOO]);
%T_MOO=[0;679.64;84.05]+(L_Moo_hat*(-15.47))-DD*L_Moo_hat;
%F1_MOO=[1;-17.6764;1.2743;-19.3393];
%F2_MOO=[0 0 0;27.1200   59.8367   38.3850;0.0000  -11.4293  130.8163;23.1300  120.6639   77.4055];

Ac_MOO=[0 1 0 0;0 0 1 0;0 0 0 1;0 -679.6 84.05 -15.47];
Bc_MOO=[0;0;0;1];
Aaa=0;
Aab=[1 0 0];
Abb=[0 1 0;0 0 1;-679.6 84.05 -15.47];
Aba=[0;0;0];
Ba=0;Bb=[0;0;1];
L_MOO=[-10 -10 -10];
Ke_MOO=acker(Abb',Aab',L_MOO)';
T_MOO=Aba-(Ke_MOO*Aaa);
D_MOO=Abb-(Ke_MOO*Aab);
R_MOO=Bb-Ke_MOO*Ba;
