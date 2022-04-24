% Solving the designing matrix in the switching function by Theorem 1

clear all
clc

% system parameters
l=2.8;
Ll=5.5;
v=-1;
bt=2;
t0=0.5;

xi=pi/2;
a=sin(xi)/xi;
E=[1 0 0 0;...
   0 1 0 0;...
   0 0 0 0;...
   0 0 0 1];
A1=[0 0    0     v*bt/(Ll*t0);...
    0 0 v*bt/t0       0      ;...
    1 0    -1    v*bt/(2*Ll) ;...
    0 0    0     -v*bt/(Ll*t0)];
A2=[0 0    0     v*bt/(Ll*t0) ;...
    0 0 v*bt/t0       0       ;...
    a 0    -1    a*v*bt/(2*Ll);...
    0 0    0     -v*bt/(Ll*t0)];
B=[0;0;0;v*bt/(l*t0)];

%  E1,Ai1,Ai2 in Theorem 1
E1=[1 0 0;...
    0 1 0;...
    0 0 0];
A11=[0 0    0   ;...
     0 0 v*bt/t0;...
     1 0   -1];
A21=[0 0    0   ;...
     0 0 v*bt/t0;...
     a 0 -1];
A12=[v*bt/(Ll*t0);0;v*bt/(2*Ll)];
A22=[v*bt/(Ll*t0);0;a*v*bt/(2*Ll)];

% U1,U2,V2 in Theorem 1
U1=[1 0 0;...
    0 1 0];
U2=[0 0 1];
V2=[0;0;1];


setlmis([]);
% Pl1
[P1,n,sP1]=lmivar(1,[2 1]);
for i=1:2
    % Pl2i
    [P2{i},n,sP2{i}]=lmivar(2,[1 2]);
    % Pl3i
    [P3{i},n,sP3{i}]=lmivar(1,[1 1]);
    % Pli
    P{i}=lmivar(3,[sP1 sP2{i}';sP2{i} sP3{i}]);
    % Phili
    Phi{i}=lmivar(2,[1 1]);
end
% Kl
K=lmivar(2,[1 2]);

lmiterm([1 1 1 P1],-1,1);

lmiterm([2 1 1 P{1}],A11,E1','s');
lmiterm([2 1 1 Phi{1}],A11*V2,U2,'s');
lmiterm([2 1 1 K],-A12,U1,'s');

lmiterm([3 1 1 P{2}],A21,E1','s');
lmiterm([3 1 1 Phi{2}],A21*V2,U2,'s');
lmiterm([3 1 1 K],-A22,U1,'s');

lmiterm([4 1 1 P{1}],A11,E1','s');
lmiterm([4 1 1 Phi{1}],A11*V2,U2,'s');
lmiterm([4 1 1 K],-A12,U1,'s');
lmiterm([4 1 1 P{2}],0.5*A11,E1','s');
lmiterm([4 1 1 Phi{2}],0.5*A11*V2,U2,'s');
lmiterm([4 1 1 K],-0.5*A12,U1,'s');
lmiterm([4 1 1 P{1}],0.5*A21,E1','s');
lmiterm([4 1 1 Phi{1}],0.5*A21*V2,U2,'s');
lmiterm([4 1 1 K],-0.5*A22,U1,'s');

lmiterm([5 1 1 P{2}],A21,E1','s');
lmiterm([5 1 1 Phi{2}],A21*V2,U2,'s');
lmiterm([5 1 1 K],-A22,U1,'s');
lmiterm([5 1 1 P{1}],0.5*A21,E1','s');
lmiterm([5 1 1 Phi{1}],0.5*A21*V2,U2,'s');
lmiterm([5 1 1 K],-0.5*A22,U1,'s');
lmiterm([5 1 1 P{2}],0.5*A11,E1','s');
lmiterm([5 1 1 Phi{2}],0.5*A11*V2,U2,'s');
lmiterm([5 1 1 K],-0.5*A12,U1,'s');

lmisys=getlmis;
[tmin xfeas]=feasp(lmisys);
if tmin<0
mP1=dec2mat(lmisys,xfeas,P1)
mK=dec2mat(lmisys,xfeas,K)
Sl=mK*inv(mP1)*U1
end



