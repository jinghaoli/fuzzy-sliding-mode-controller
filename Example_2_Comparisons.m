% Comparison among classical (Theorem 1), PDC fuzzy (Corollary 1) and Non-PDC (Theorem 3) fuzzy-based methods
%        it shows that the nonpdc is better
%------------------------------------------------------------------------
clear all
clc

figure;

for beta=-3:0.2:0
    for gamma=0:0.2:3
        E=[1 0 0;0 1 0;0 0 0];
        A1{1}=[-1 beta 1;-0.5 -0.5 -1;0 0 -1];
        A2{1}=[-0.5;1;1];
        A1{2}=[gamma 0.5 0;0.5 -1 0;1 0 -1];
        A2{2}=[1;0.5;0.5];
        A1{3}=[-1 beta 0;gamma -1 1;0 1 -1];
        A2{3}=[1;0;1];

        U1=[1 0 0;0 1 0];
        U2=[0 0 1];
        V2=[0;0;1];
        
        
        %---------------------------------------------------------------------------
        %      Classical linear switching function-based method
        %---------------------------------------------------------------------------
        setlmis([]);
        
        %Pc1
        [Pc1,n,sPc1]=lmivar(1,[2 1]);
        for i=1:3
            %Pc2i
            [Pc2{i},n,sPc2{i}]=lmivar(2,[1 2]);
            %Pc3i
            [Pc3{i},n,sPc3{i}]=lmivar(1,[1 1]);
            %Pci
            Pc{i}=lmivar(3,[sPc1 sPc2{i}';sPc2{i} sPc3{i}]);
            %Phici
            Phic{i}=lmivar(2,[1 1]);
        end
        %Kc
        Kc=lmivar(2,[1 2]);
        
        lmiterm([1 1 1 Pc{1}],A1{1},E','s');
        lmiterm([1 1 1 Phic{1}],A1{1}*V,U,'s');
        lmiterm([1 1 1 Kc],-A2{1},U1,'s');
        
        lmiterm([2 1 1 Pc{2}],A1{2},E','s');
        lmiterm([2 1 1 Phic{2}],A1{2}*V,U,'s');
        lmiterm([2 1 1 Kc],-A2{2},U1,'s');
        
        lmiterm([3 1 1 Pc{3}],A1{3},E','s');
        lmiterm([3 1 1 Phic{3}],A1{3}*V,U,'s');
        lmiterm([3 1 1 Kc],-A2{3},U1,'s');
        
        lmiterm([4 1 1 Pc{1}],0.5*A1{1},E','s');
        lmiterm([4 1 1 Phic{1}],0.5*A1{1}*V,U,'s');
        lmiterm([4 1 1 Kc],-A2{1},0.5*U1,'s');
        lmiterm([4 1 1 Pc{2}],0.5*A1{1},E','s');
        lmiterm([4 1 1 Phic{2}],0.5*A1{1}*V,U,'s');
        lmiterm([4 1 1 Kc],-0.5*A2{1},U1,'s');
        lmiterm([4 1 1 Pc{1}],0.5*A1{2},E','s');
        lmiterm([4 1 1 Phic{1}],0.5*A1{2}*V,U,'s');
        lmiterm([4 1 1 Kc],-0.5*A2{2},U1,'s');
        
        lmiterm([5 1 1 Pc{1}],0.5*A1{1},E','s');
        lmiterm([5 1 1 Phic{1}],0.5*A1{1}*V,U,'s');
        lmiterm([5 1 1 Kc],-A2{1},0.5*U1,'s');
        lmiterm([5 1 1 Pc{3}],0.5*A1{1},E','s');
        lmiterm([5 1 1 Phic{3}],0.5*A1{1}*V,U,'s');
        lmiterm([5 1 1 Kc],-0.5*A2{1},U1,'s');
        lmiterm([5 1 1 Pc{1}],0.5*A1{3},E','s');
        lmiterm([5 1 1 Phic{1}],0.5*A1{3}*V,U,'s');
        lmiterm([5 1 1 Kc],-0.5*A2{3},U1,'s');
        
        lmiterm([6 1 1 Pc{2}],0.5*A1{2},E','s');
        lmiterm([6 1 1 Phic{2}],0.5*A1{2}*V,U,'s');
        lmiterm([6 1 1 Kc],-A2{2},0.5*U1,'s');
        lmiterm([6 1 1 Pc{1}],0.5*A1{2},E','s');
        lmiterm([6 1 1 Phic{1}],0.5*A1{2}*V,U,'s');
        lmiterm([6 1 1 Kc],-0.5*A2{2},U1,'s');
        lmiterm([6 1 1 Pc{2}],0.5*A1{1},E','s');
        lmiterm([6 1 1 Phic{2}],0.5*A1{1}*V,U,'s');
        lmiterm([6 1 1 Kc],-0.5*A2{1},U1,'s');
        
        lmiterm([7 1 1 Pc{2}],0.5*A1{2},E','s');
        lmiterm([7 1 1 Phic{2}],0.5*A1{2}*V,U,'s');
        lmiterm([7 1 1 Kc],-A2{2},0.5*U1,'s');
        lmiterm([7 1 1 Pc{3}],0.5*A1{2},E','s');
        lmiterm([7 1 1 Phic{3}],0.5*A1{2}*V,U,'s');
        lmiterm([7 1 1 Kc],-0.5*A2{2},U1,'s');
        lmiterm([7 1 1 Pc{2}],0.5*A1{3},E','s');
        lmiterm([7 1 1 Phic{2}],0.5*A1{3}*V,U,'s');
        lmiterm([7 1 1 Kc],-0.5*A2{3},U1,'s');
        
        lmiterm([8 1 1 Pc{3}],0.5*A1{3},E','s');
        lmiterm([8 1 1 Phic{3}],0.5*A1{3}*V,U,'s');
        lmiterm([8 1 1 Kc],-0.5*A2{3},U1,'s');
        lmiterm([8 1 1 Pc{1}],0.5*A1{3},E','s');
        lmiterm([8 1 1 Phic{1}],0.5*A1{3}*V,U,'s');
        lmiterm([8 1 1 Kc],-0.5*A2{3},U1,'s');
        lmiterm([8 1 1 Pc{3}],0.5*A1{1},E','s');
        lmiterm([8 1 1 Phic{3}],0.5*A1{1}*V,U,'s');
        lmiterm([8 1 1 Kc],-0.5*A2{1},U1,'s');
        
        lmiterm([9 1 1 Pc{3}],0.5*A1{3},E','s');
        lmiterm([9 1 1 Phic{3}],0.5*A1{3}*V,U,'s');
        lmiterm([9 1 1 Kc],-A2{3},0.5*U1,'s');
        lmiterm([9 1 1 Pc{2}],0.5*A1{3},E','s');
        lmiterm([9 1 1 Phic{2}],0.5*A1{3}*V,U,'s');
        lmiterm([9 1 1 Kc],-0.5*A2{3},U1,'s');
        lmiterm([9 1 1 Pc{3}],0.5*A1{2},E','s');
        lmiterm([9 1 1 Phic{3}],0.5*A1{2}*V,U,'s');
        lmiterm([9 1 1 Kc],-0.5*A2{2},U1,'s');
        
        lmiterm([10 1 1 Pc1],-1,1);
        
        lmisys=getlmis;
        [tminc xfeas]=feasp(lmisys);
        if tminc<0
            plot(beta,gamma,'bo');hold on
        end
        
        %----------------------------------------------------------------------------
        %   PDC fuzzy switching function-based method
        %---------------------------------------------------------------------------
        setlmis([]);
        
        %Pp1
        [Pp1,n,sPp1]=lmivar(1,[2 1]);
        for i=1:3
            %Pp2i
            [Pp2{i},n,sPp2{i}]=lmivar(2,[1 2]);
            %Pp3i
            [Pp3{i},n,sPp3{i}]=lmivar(1,[1 1]);
            %Ppi
            Pp{i}=lmivar(3,[sPp1 sPp2{i}';sPp2{i} sPp3{i}]);
            %Phipi
            Phip{i}=lmivar(2,[1 1]);
            %Kpi
            Kp{i}=lmivar(2,[1 3]);
        end
        
        lmiterm([1 1 1 Pp{1}],A1{1},E','s');
        lmiterm([1 1 1 Phip{1}],A1{1}*V,U,'s');
        lmiterm([1 1 1 Kp{1}],-A2{1},E','s');
        
        lmiterm([2 1 1 Pp{2}],A1{2},E','s');
        lmiterm([2 1 1 Phip{2}],A1{2}*V,U,'s');
        lmiterm([2 1 1 Kp{2}],-A2{2},E','s');
        
        lmiterm([3 1 1 Pp{3}],A1{3},E','s');
        lmiterm([3 1 1 Phip{3}],A1{3}*V,U,'s');
        lmiterm([3 1 1 Kp{3}],-A2{3},E','s');
        
        lmiterm([4 1 1 Pp{1}],0.5*A1{1},E','s');
        lmiterm([4 1 1 Phip{1}],0.5*A1{1}*V,U,'s');
        lmiterm([4 1 1 Kp{1}],-0.5*A2{1},E','s');
        lmiterm([4 1 1 Pp{2}],0.5*A1{1},E','s');
        lmiterm([4 1 1 Phip{2}],0.5*A1{1}*V,U,'s');
        lmiterm([4 1 1 Kp{2}],-0.5*A2{1},E','s');
        lmiterm([4 1 1 Pp{1}],0.5*A1{2},E','s');
        lmiterm([4 1 1 Phip{1}],0.5*A1{2}*V,U,'s');
        lmiterm([4 1 1 Kp{1}],-0.5*A2{2},E','s');
        
        lmiterm([5 1 1 Pp{1}],0.5*A1{1},E','s');
        lmiterm([5 1 1 Phip{1}],0.5*A1{1}*V,U,'s');
        lmiterm([5 1 1 Kp{1}],-0.5*A2{1},E','s');
        lmiterm([5 1 1 Pp{3}],0.5*A1{1},E','s');
        lmiterm([5 1 1 Phip{3}],0.5*A1{1}*V,U,'s');
        lmiterm([5 1 1 Kp{3}],-0.5*A2{1},E','s');
        lmiterm([5 1 1 Pp{1}],0.5*A1{3},E','s');
        lmiterm([5 1 1 Phip{1}],0.5*A1{3}*V,U,'s');
        lmiterm([5 1 1 Kp{1}],-0.5*A2{3},E','s');
        
        lmiterm([6 1 1 Pp{2}],0.5*A1{2},E','s');
        lmiterm([6 1 1 Phip{2}],0.5*A1{2}*V,U,'s');
        lmiterm([6 1 1 Kp{2}],-0.5*A2{2},E','s');
        lmiterm([6 1 1 Pp{1}],0.5*A1{2},E','s');
        lmiterm([6 1 1 Phip{1}],0.5*A1{2}*V,U,'s');
        lmiterm([6 1 1 Kp{1}],-0.5*A2{2},E','s');
        lmiterm([6 1 1 Pp{2}],0.5*A1{1},E','s');
        lmiterm([6 1 1 Phip{2}],0.5*A1{1}*V,U,'s');
        lmiterm([6 1 1 Kp{2}],-0.5*A2{1},E','s');
        
        lmiterm([7 1 1 Pp{2}],0.5*A1{2},E','s');
        lmiterm([7 1 1 Phip{2}],0.5*A1{2}*V,U,'s');
        lmiterm([7 1 1 Kp{2}],-0.5*A2{2},E','s');
        lmiterm([7 1 1 Pp{3}],0.5*A1{2},E','s');
        lmiterm([7 1 1 Phip{3}],0.5*A1{2}*V,U,'s');
        lmiterm([7 1 1 Kp{3}],-0.5*A2{2},E','s');
        lmiterm([7 1 1 Pp{2}],0.5*A1{3},E','s');
        lmiterm([7 1 1 Phip{2}],0.5*A1{3}*V,U,'s');
        lmiterm([7 1 1 Kp{2}],-0.5*A2{3},E','s');
        
        lmiterm([8 1 1 Pp{3}],0.5*A1{3},E','s');
        lmiterm([8 1 1 Phip{3}],0.5*A1{3}*V,U,'s');
        lmiterm([8 1 1 Kp{3}],-0.5*A2{3},E','s');
        lmiterm([8 1 1 Pp{1}],0.5*A1{3},E','s');
        lmiterm([8 1 1 Phip{1}],0.5*A1{3}*V,U,'s');
        lmiterm([8 1 1 Kp{1}],-0.5*A2{3},E','s');
        lmiterm([8 1 1 Pp{3}],0.5*A1{1},E','s');
        lmiterm([8 1 1 Phip{3}],0.5*A1{1}*V,U,'s');
        lmiterm([8 1 1 Kp{3}],-0.5*A2{1},E','s');
        
        lmiterm([9 1 1 Pp{3}],0.5*A1{3},E','s');
        lmiterm([9 1 1 Phip{3}],0.5*A1{3}*V,U,'s');
        lmiterm([9 1 1 Kp{3}],-0.5*A2{3},E','s');
        lmiterm([9 1 1 Pp{2}],0.5*A1{3},E','s');
        lmiterm([9 1 1 Phip{2}],0.5*A1{3}*V,U,'s');
        lmiterm([9 1 1 Kp{2}],-0.5*A2{3},E','s');
        lmiterm([9 1 1 Pp{3}],0.5*A1{2},E','s');
        lmiterm([9 1 1 Phip{3}],0.5*A1{2}*V,U,'s');
        lmiterm([9 1 1 Kp{3}],-0.5*A2{2},E','s');
        
        lmiterm([10 1 1 Pp1],-1,1);
        
        lmisys=getlmis;
        [tminp xfeas]=feasp(lmisys);
        if tminp<0
            plot(beta,gamma,'r*');hold on
        end
        
        %----------------------------------------------------------------------------
        %   Non-PDC fuzzy switching function-based method
        %   Solved by YALMIP and SeDuMi due to the non-strict matrix
        %   inequalities
        %---------------------------------------------------------------------------
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
            plot(beta,gamma,'k+');hold on
        end
    end
end
xlabel('a','FontSize',25);
ylabel('b','FontSize',25);
set(gca,'FontSize',25);




