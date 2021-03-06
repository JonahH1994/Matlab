function A = AMatrix(I1,I2,I3,d1,d2,d3,l1,l2,l3,m1,m2,m3,th1,th2,th3,x1,x2,x3,y1,y2,y3)
%AMATRIX
%    A = AMATRIX(I1,I2,I3,D1,D2,D3,L1,L2,L3,M1,M2,M3,TH1,TH2,TH3,X1,X2,X3,Y1,Y2,Y3)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    12-Dec-2016 10:56:45

t2 = 1.0./l1;
t3 = d1.*m1.*t2;
t4 = 1.0./l2;
t5 = d2.*t4;
t6 = t5-1.0;
t7 = d2.*m2.*t4;
t8 = 1.0./l3;
t9 = d3.*t8;
t10 = t9-1.0;
t11 = d3.*m3.*t8;
t12 = d1.*t2.*x1;
t13 = y1-y2;
t14 = x1-x2;
t15 = d2.*t4.*t14;
t16 = cos(th1);
t17 = sin(th1);
t18 = t16.^2;
t19 = t17.^2;
t20 = t18+t19;
t21 = 1.0./t20;
t22 = cos(th2);
t23 = sin(th2);
t24 = t22.^2;
t25 = t23.^2;
t26 = t24+t25;
t27 = 1.0./t26;
t28 = cos(th3);
t29 = sin(th3);
t30 = l1.*t16.*t21;
t31 = t28.^2;
t32 = t29.^2;
t33 = t31+t32;
t34 = 1.0./t33;
A = reshape([t3,0.0,-m2.*t6,0.0,0.0,0.0,0.0,0.0,0.0,-d1.*t2,0.0,t6,0.0,0.0,0.0,0.0,0.0,t7,0.0,-m3.*t10,0.0,0.0,0.0,0.0,0.0,0.0,-t5,0.0,t10,0.0,0.0,0.0,0.0,0.0,t11,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t9,0.0,0.0,t3,0.0,-m2.*t6,0.0,0.0,0.0,0.0,0.0,0.0,-d1.*t2,0.0,t6,0.0,0.0,0.0,0.0,0.0,t7,0.0,-m3.*t10,0.0,0.0,0.0,0.0,0.0,0.0,-t5,0.0,t10,0.0,0.0,0.0,0.0,0.0,t11,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t9,0.0,0.0,0.0,0.0,0.0,0.0,I1,0.0,0.0,-d1.*t17.*t21,d1.*t16.*t21,-l1.*t17.*t21,t30,-l1.*t17.*t21,t30,0.0,0.0,0.0,0.0,0.0,0.0,0.0,I2,0.0,0.0,0.0,-d2.*t23.*t27,d2.*t22.*t27,-l2.*t23.*t27,l2.*t22.*t27,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,I3,0.0,0.0,0.0,0.0,-d3.*t29.*t34,d3.*t28.*t34,-1.0,0.0,0.0,0.0,0.0,0.0,-d1.*t2.*y1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,t12,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,y1-d1.*t2.*y1,-d2.*t4.*t13,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,t12-x1,t15,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,-y1+y2+d2.*t4.*t13,-d3.*t8.*(y2-y3),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,-t15+x1-x2,d3.*t8.*(x2-x3),0.0,0.0,0.0,0.0,0.0,0.0],[15,15]);
