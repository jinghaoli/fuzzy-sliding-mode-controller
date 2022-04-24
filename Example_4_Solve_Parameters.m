% Solve the design paramters in switching function by Corollary 1

clear all
clc

% System simulation parameters
Jm=0.126*10^-3;
alpha=0.001;
R=0.1;
L=92;
ke=0.076;
l=0.2;
m=0.1;
km=0.0446;
g=9.81;

% Paramters in the T-S fuzzy systems
xi=pi/3;
a1=1;
a2=sin(xi)/xi;
disp('**********************************************************');
disp('*********Paramters in the T-S fuzzy systems***************');
A1=[0 1 0;-m*g*l*a1/(Jm+m*l*l) -alpha/(Jm+m*l*l) km/(Jm+m*l*l);...
    0 -ke/L -R/L];
A2=[0 1 0;-m*g*l*a2/(Jm+m*l*l) -alpha/(Jm+m*l*l) km/(Jm+m*l*l);...
    0 -ke/L -R/L];
B=[0;0;1/L];
A11=[0 1;-m*g*l*a1/(Jm+m*l*l) -alpha/(Jm+m*l*l)]
A21=[0 1;-m*g*l*a2/(Jm+m*l*l) -alpha/(Jm+m*l*l)]
A12=[0;km/(Jm+m*l*l)]
A22=[0;km/(Jm+m*l*l)]
A13=[0 -ke/L]
A23=[0 -ke/L]
A14=-R/L
A24=-R/L
B1=1/L
disp('**********************************************************');

bA11=A11+3*[1 0;0 1];
bA21=A21+3*[1 0;0 1];

% Solve the conditions in Corollary 1
setlmis([]);
%Pp
Pp=lmivar(1,[2 1]);
Kp1=lmivar(2,[1 2]);
Kp2=lmivar(2,[1 2]);

lmiterm([1 1 1 Pp],bA11,1,'s');
lmiterm([1 1 1 Kp1],-A12,1,'s');

lmiterm([2 1 1 Pp],bA21,1,'s');
lmiterm([2 1 1 Kp2],-A22,1,'s');

lmiterm([3 1 1 Pp],bA11,1,'s');
lmiterm([3 1 1 Kp1],-A12,1,'s');
lmiterm([3 1 1 Pp],0.5*bA11,1,'s');
lmiterm([3 1 1 Kp2],-0.5*A12,1,'s');
lmiterm([3 1 1 Pp],0.5*bA21,1,'s');
lmiterm([3 1 1 Kp1],-0.5*A22,1,'s');

lmiterm([4 1 1 Pp],bA21,1,'s');
lmiterm([4 1 1 Kp2],-A22,1,'s');
lmiterm([4 1 1 Pp],0.5*bA21,1,'s');
lmiterm([4 1 1 Kp1],-0.5*A22,1,'s');
lmiterm([4 1 1 Pp],0.5*bA11,1,'s');
lmiterm([4 1 1 Kp2],-0.5*A12,1,'s');

lmiterm([5 1 1 Pp],-1,1);

lmisys=getlmis;
[tminc xfeas]=feasp(lmisys);
if tminc<0
mPp=dec2mat(lmisys,xfeas,Pp);
mKp1=dec2mat(lmisys,xfeas,Kp1);
mKp2=dec2mat(lmisys,xfeas,Kp2);
disp('**** The design parameters in the switching function ****');
Sp1=mKp1*inv(mPp)
Sp2=mKp2*inv(mPp)
end

