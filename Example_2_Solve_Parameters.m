%----------------------------------------------------------------------------
%   Non-PDC fuzzy switching function-based method
%   Solved by YALMIP and SeDuMi due to the non-strict matrix
%   inequalities
%---------------------------------------------------------------------------
clear all
clc

beta=-2;
gamma=1;
E=[1 0 0;0 1 0;0 0 0];
A1{1}=[-1 beta 1;-0.5 -0.5 -1;0 0 -1];
A2{1}=[-0.5;1;1];
A1{2}=[gamma 0.5 0;0.5 -1 0;1 0 -1];
A2{2}=[1;0.5;0.5];
A1{3}=[-1 beta 0;gamma -1 1;0 1 -1];
A2{3}=[1;0;1];

U1=[1 0 0;0 1 0];
U=[0 0 1];
V=[0;0;1];
phi=-2;

bl=10^-12; % the constant is used to ensure that the strict matrix inequalities hold
for i=1:3
    %Pn1
    Pn1{i}=sdpvar(2,2);
    %Pn2
    Pn2{i}=sdpvar(1,2,'full');
    % Pn3 is not defined since it will not appear in the matrix
    % inequalities
    %Phin
    Phin{i}=sdpvar(1,1);
    %Kni
    Kn{i}=sdpvar(1,3,'full');
end
Xn=sdpvar(2,2);

Fn=[A1{1}*[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]+[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]'*A1{1}'-A2{1}*Kn{1}*E'-E*Kn{1}'*A2{1}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'<=-bl*eye(3)];
Fn=[Fn,A1{2}*[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]+[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]'*A1{2}'-A2{2}*Kn{2}*E'-E*Kn{2}'*A2{2}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'<=-bl*eye(3)];
Fn=[Fn,A1{3}*[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]+[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]'*A1{3}'-A2{3}*Kn{3}*E'-E*Kn{3}'*A2{3}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'<=-bl*eye(3)];
Fn=[Fn,0.5*(A1{1}*[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]+[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]'*A1{1}'-A2{1}*Kn{1}*E'-E*Kn{1}'*A2{1}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'+A1{1}*[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]+[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]'*A1{1}'-A2{1}*Kn{2}*E'-E*Kn{2}'*A2{1}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'+A1{2}*[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]+[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]'*A1{2}'-A2{2}*Kn{1}*E'-E*Kn{1}'*A2{2}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E')<=-bl*eye(3)];
Fn=[Fn,0.5*(A1{1}*[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]+[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]'*A1{1}'-A2{1}*Kn{1}*E'-E*Kn{1}'*A2{1}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'+A1{1}*[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]+[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]'*A1{1}'-A2{1}*Kn{3}*E'-E*Kn{3}'*A2{1}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'+A1{3}*[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]+[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]'*A1{3}'-A2{3}*Kn{1}*E'-E*Kn{1}'*A2{3}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E')<=-bl*eye(3)];
Fn=[Fn,0.5*(A1{2}*[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]+[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]'*A1{2}'-A2{2}*Kn{2}*E'-E*Kn{2}'*A2{2}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'+A1{2}*[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]+[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]'*A1{2}'-A2{2}*Kn{1}*E'-E*Kn{1}'*A2{2}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'+A1{1}*[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]+[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]'*A1{1}'-A2{1}*Kn{2}*E'-E*Kn{2}'*A2{1}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E')<=-bl*eye(3)];
Fn=[Fn,0.5*(A1{2}*[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]+[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]'*A1{2}'-A2{2}*Kn{2}*E'-E*Kn{2}'*A2{2}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'+A1{2}*[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]+[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]'*A1{2}'-A2{2}*Kn{3}*E'-E*Kn{3}'*A2{2}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'+A1{3}*[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]+[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]'*A1{3}'-A2{3}*Kn{2}*E'-E*Kn{2}'*A2{3}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E')<=-bl*eye(3)];
Fn=[Fn,0.5*(A1{3}*[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]+[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]'*A1{3}'-A2{3}*Kn{3}*E'-E*Kn{3}'*A2{3}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'+A1{3}*[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]+[Pn1{1} zeros(2,1);Pn2{1} Phin{1}]'*A1{3}'-A2{3}*Kn{1}*E'-E*Kn{1}'*A2{3}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'+A1{1}*[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]+[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]'*A1{1}'-A2{1}*Kn{3}*E'-E*Kn{3}'*A2{1}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E')<=-bl*eye(3)];
Fn=[Fn,0.5*(A1{3}*[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]+[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]'*A1{3}'-A2{3}*Kn{3}*E'-E*Kn{3}'*A2{3}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'+A1{3}*[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]+[Pn1{2} zeros(2,1);Pn2{2} Phin{2}]'*A1{3}'-A2{3}*Kn{2}*E'-E*Kn{2}'*A2{3}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E'+A1{2}*[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]+[Pn1{3} zeros(2,1);Pn2{3} Phin{3}]'*A1{2}'-A2{2}*Kn{3}*E'-E*Kn{3}'*A2{2}'-phi*E*[Pn1{1}+Pn1{2}+Pn1{3}+3*Xn zeros(2,1);zeros(1,2) 0]*E')<=-bl*eye(3)];
Fn=[Fn,Pn1{1}>=bl*eye(2)];
Fn=[Fn,Pn1{2}>=bl*eye(2)];
Fn=[Fn,Pn1{3}>=bl*eye(2)];
Fn=[Fn,Pn1{1}+Xn>=0];
Fn=[Fn,Pn1{2}+Xn>=0];
Fn=[Fn,Pn1{3}+Xn>=0];

diagnosticsn=optimize(Fn);

if diagnosticsn.problem==0
    mPn11=value(Pn1{1});
    mPn12=value(Pn1{2});
    mPn13=value(Pn1{3});
    mPn21=value(Pn2{1});
    mPn22=value(Pn2{2});
    mPn23=value(Pn2{3});
    mPhin1=value(Phin{1});
    mPhin2=value(Phin{2});
    mPhin3=value(Phin{3});
    mYn1=[mPn11 zeros(2,1);mPn21 mPhin1]
    mYn2=[mPn12 zeros(2,1);mPn22 mPhin2]
    mYn3=[mPn13 zeros(2,1);mPn23 mPhin3]
    mKn1=value(Kn{1})
    mKn2=value(Kn{2})
    mKn3=value(Kn{3})
end


