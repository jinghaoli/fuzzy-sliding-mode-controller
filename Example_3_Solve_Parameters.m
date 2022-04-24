% Solving the designing parameter in the switching function by Theorem 1
clear all
clc

% System parameters
m=2;
M=8;
l=0.5;
k=0.5;
xi=pi/3;
a=cos(xi);
A11=[0    0     0;...
     0    0     1;...
     0 -k/(m+M) 0];
A21=[0    0     0;...
     0    0     1;...
     0 -k/(m+M) 0];
A12=[1;-m*l/(m+M);0];
A22=[1;-a*m*l/(m+M);0];

% Let the real part of all the eigenvalues of the closed-loop system is less
% than 0.1
bA11=A11+0.1*[1 0 0;0 1 0;0 0 1];
bA21=A21+0.1*[1 0 0;0 1 0;0 0 1];


setlmis([]);
%Pc
Pc=lmivar(1,[3 1]);
%Kc
Kc=lmivar(2,[1 3]);

lmiterm([1 1 1 Pc],bA11,1,'s');
lmiterm([1 1 1 Kc],-A12,1,'s');

lmiterm([2 1 1 Pc],bA21,1,'s');
lmiterm([2 1 1 Kc],-A22,1,'s');

lmiterm([3 1 1 Pc],-1,1);
lmiterm([3 1 1 0],0);

lmisys=getlmis;
[tminc xfeas]=feasp(lmisys);
if tminc<0
mPc=dec2mat(lmisys,xfeas,Pc)
mKc=dec2mat(lmisys,xfeas,Kc)
mSc=mKc*inv(mPc)
end

