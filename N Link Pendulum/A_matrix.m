function A = A_matrix(in1,in2,in3,in4)
%A_MATRIX
%    A = A_MATRIX(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    16-Dec-2016 13:05:55

I1 = in1(:,1);
I2 = in1(:,2);
I3 = in1(:,3);
I4 = in1(:,4);
I5 = in1(:,5);
I6 = in1(:,6);
I7 = in1(:,7);
I8 = in1(:,8);
I9 = in1(:,9);
l1 = in2(:,1);
l2 = in2(:,2);
l3 = in2(:,3);
l4 = in2(:,4);
l5 = in2(:,5);
l6 = in2(:,6);
l7 = in2(:,7);
l8 = in2(:,8);
l9 = in2(:,9);
m1 = in3(:,1);
m2 = in3(:,2);
m3 = in3(:,3);
m4 = in3(:,4);
m5 = in3(:,5);
m6 = in3(:,6);
m7 = in3(:,7);
m8 = in3(:,8);
m9 = in3(:,9);
th1 = in4(:,1);
th2 = in4(:,2);
th3 = in4(:,3);
th4 = in4(:,4);
th5 = in4(:,5);
th6 = in4(:,6);
th7 = in4(:,7);
th8 = in4(:,8);
th9 = in4(:,9);
t2 = sin(th1);
t3 = cos(th1);
t4 = l1.*t2;
t5 = sin(th2);
t6 = l2.*t5;
t7 = sin(th3);
t8 = l1.*t3;
t9 = cos(th2);
t10 = l2.*t9;
t11 = cos(th3);
t12 = l3.*t7;
t13 = sin(th4);
t14 = l4.*t13;
t15 = sin(th5);
t16 = l5.*t15;
t17 = sin(th6);
t18 = l3.*t11;
t19 = cos(th4);
t20 = l4.*t19;
t21 = cos(th5);
t22 = l5.*t21;
t23 = cos(th6);
t24 = l6.*t23;
t25 = cos(th7);
t26 = l6.*t17;
t27 = sin(th7);
t28 = l7.*t25;
t29 = cos(th8);
t30 = l7.*t27;
t31 = sin(th8);
t32 = l3.*t7.*(1.0./2.0);
t33 = t4+t6+t32;
t34 = l3.*t11.*(1.0./2.0);
t35 = t8+t10+t34;
t36 = l6.*t17.*(1.0./2.0);
t37 = t4+t6+t12+t14+t16+t36;
t38 = l6.*t23.*(1.0./2.0);
t39 = t8+t10+t18+t20+t22+t38;
t40 = l8.*t31;
t41 = sin(th9);
t42 = l9.*t41.*(1.0./2.0);
t43 = t4+t6+t12+t14+t16+t26+t30+t40+t42;
t44 = l8.*t29;
t45 = cos(th9);
t46 = l9.*t45.*(1.0./2.0);
t47 = t8+t10+t18+t20+t22+t24+t28+t44+t46;
t48 = l4.*t19.*(1.0./2.0);
t49 = t8+t10+t18+t48;
t50 = l4.*t13.*(1.0./2.0);
t51 = t4+t6+t12+t50;
t52 = l7.*t25.*(1.0./2.0);
t53 = t8+t10+t18+t20+t22+t24+t52;
t54 = l7.*t27.*(1.0./2.0);
t55 = t4+t6+t12+t14+t16+t26+t54;
t56 = l2.*t9.*(1.0./2.0);
t57 = t8+t56;
t58 = l2.*t5.*(1.0./2.0);
t59 = t4+t58;
t60 = l5.*t21.*(1.0./2.0);
t61 = t8+t10+t18+t20+t60;
t62 = l5.*t15.*(1.0./2.0);
t63 = t4+t6+t12+t14+t62;
t64 = l8.*t29.*(1.0./2.0);
t65 = t8+t10+t18+t20+t22+t24+t28+t64;
t66 = l8.*t31.*(1.0./2.0);
t67 = t4+t6+t12+t14+t16+t26+t30+t66;
t68 = t6+t12+t50;
t69 = t10+t18+t48;
t70 = t6+t12+t14+t16+t26+t54;
t71 = t10+t18+t20+t22+t24+t52;
t72 = t10+t18+t20+t60;
t73 = t6+t12+t14+t62;
t74 = t10+t18+t20+t22+t24+t28+t64;
t75 = t6+t12+t14+t16+t26+t30+t66;
t76 = t10+t34;
t77 = t6+t32;
t78 = t10+t18+t20+t22+t38;
t79 = t6+t12+t14+t16+t36;
t80 = t10+t18+t20+t22+t24+t28+t44+t46;
t81 = t6+t12+t14+t16+t26+t30+t40+t42;
t82 = t12+t14+t62;
t83 = t18+t20+t60;
t84 = t12+t14+t16+t26+t30+t66;
t85 = t18+t20+t22+t24+t28+t64;
t86 = t18+t20+t22+t38;
t87 = t12+t14+t16+t36;
t88 = t18+t20+t22+t24+t28+t44+t46;
t89 = t12+t14+t16+t26+t30+t40+t42;
t90 = t18+t48;
t91 = t12+t50;
t92 = t18+t20+t22+t24+t52;
t93 = t12+t14+t16+t26+t54;
t94 = l3.^2;
t95 = t14+t16+t36;
t96 = t20+t22+t38;
t97 = t14+t16+t26+t30+t40+t42;
t98 = t20+t22+t24+t28+t44+t46;
t99 = t20+t22+t24+t52;
t100 = t14+t16+t26+t54;
t101 = t20+t60;
t102 = t14+t62;
t103 = t20+t22+t24+t28+t64;
t104 = t14+t16+t26+t30+t66;
t105 = l4.^2;
t106 = t16+t26+t54;
t107 = t22+t24+t52;
t108 = t22+t24+t28+t64;
t109 = t16+t26+t30+t66;
t110 = t22+t38;
t111 = t16+t36;
t112 = t22+t24+t28+t44+t46;
t113 = t16+t26+t30+t40+t42;
t114 = l5.^2;
t115 = t26+t30+t66;
t116 = t24+t28+t64;
t117 = t24+t28+t44+t46;
t118 = t26+t30+t40+t42;
t119 = t24+t52;
t120 = t26+t54;
t121 = l6.^2;
t122 = t30+t40+t42;
t123 = t28+t44+t46;
t124 = t28+t64;
t125 = t30+t66;
t126 = l7.^2;
t127 = t44+t46;
t128 = t40+t42;
t129 = l8.^2;
t130 = l9.^2;
A = reshape([I1+m3.*(l1.*t2.*t33+l1.*t3.*t35)+m6.*(l1.*t2.*t37+l1.*t3.*t39)+m9.*(l1.*t2.*t43+l1.*t3.*t47)+m4.*(l1.*t3.*t49+l1.*t2.*t51)+m7.*(l1.*t3.*t53+l1.*t2.*t55)+m2.*(l1.*t3.*t57+l1.*t2.*t59)+m5.*(l1.*t3.*t61+l1.*t2.*t63)+m8.*(l1.*t3.*t65+l1.*t2.*t67)+l1.^2.*m1.*(1.0./4.0),m4.*(l1.*t2.*t68+l1.*t3.*t69)+m7.*(l1.*t2.*t70+l1.*t3.*t71)+m5.*(l1.*t2.*t73+l1.*t3.*t72)+m3.*(l1.*t2.*t77+l1.*t3.*t76)+m8.*(l1.*t2.*t75+l1.*t3.*t74)+m6.*(l1.*t2.*t79+l1.*t3.*t78)+m9.*(l1.*t2.*t81+l1.*t3.*t80)+l1.*l2.*m2.*cos(th1-th2).*(1.0./2.0),m3.*(l1.*l3.*t2.*t7.*(1.0./2.0)+l1.*l3.*t3.*t11.*(1.0./2.0))+m5.*(l1.*t2.*t82+l1.*t3.*t83)+m8.*(l1.*t2.*t84+l1.*t3.*t85)+m6.*(l1.*t2.*t87+l1.*t3.*t86)+m4.*(l1.*t2.*t91+l1.*t3.*t90)+m9.*(l1.*t2.*t89+l1.*t3.*t88)+m7.*(l1.*t2.*t93+l1.*t3.*t92),m4.*(l1.*l4.*t2.*t13.*(1.0./2.0)+l1.*l4.*t3.*t19.*(1.0./2.0))+m6.*(l1.*t2.*t95+l1.*t3.*t96)+m9.*(l1.*t2.*t97+l1.*t3.*t98)+m7.*(l1.*t2.*t100+l1.*t3.*t99)+m5.*(l1.*t2.*t102+l1.*t3.*t101)+m8.*(l1.*t2.*t104+l1.*t3.*t103),m5.*(l1.*l5.*t2.*t15.*(1.0./2.0)+l1.*l5.*t3.*t21.*(1.0./2.0))+m7.*(l1.*t2.*t106+l1.*t3.*t107)+m8.*(l1.*t2.*t109+l1.*t3.*t108)+m6.*(l1.*t2.*t111+l1.*t3.*t110)+m9.*(l1.*t2.*t113+l1.*t3.*t112),m6.*(l1.*l6.*t2.*t17.*(1.0./2.0)+l1.*l6.*t3.*t23.*(1.0./2.0))+m8.*(l1.*t2.*t115+l1.*t3.*t116)+m9.*(l1.*t2.*t118+l1.*t3.*t117)+m7.*(l1.*t2.*t120+l1.*t3.*t119),m7.*(l1.*l7.*t3.*t25.*(1.0./2.0)+l1.*l7.*t2.*t27.*(1.0./2.0))+m9.*(l1.*t2.*t122+l1.*t3.*t123)+m8.*(l1.*t2.*t125+l1.*t3.*t124),m8.*(l1.*l8.*t3.*t29.*(1.0./2.0)+l1.*l8.*t2.*t31.*(1.0./2.0))+m9.*(l1.*t2.*t128+l1.*t3.*t127),l1.*l9.*t2.*t41.*(1.0./2.0)+l1.*l9.*t3.*t45.*(1.0./2.0),I2+m3.*(l2.*t5.*t33+l2.*t9.*t35)+m6.*(l2.*t5.*t37+l2.*t9.*t39)+m9.*(l2.*t5.*t43+l2.*t9.*t47)+m4.*(l2.*t5.*t51+l2.*t9.*t49)+m7.*(l2.*t5.*t55+l2.*t9.*t53)+m2.*(l2.*t5.*t59.*(1.0./2.0)+l2.*t9.*t57.*(1.0./2.0))+m5.*(l2.*t5.*t63+l2.*t9.*t61)+m8.*(l2.*t5.*t67+l2.*t9.*t65),I2+m4.*(l2.*t5.*t68+l2.*t9.*t69)+m7.*(l2.*t5.*t70+l2.*t9.*t71)+m5.*(l2.*t5.*t73+l2.*t9.*t72)+m3.*(l2.*t5.*t77+l2.*t9.*t76)+m8.*(l2.*t5.*t75+l2.*t9.*t74)+m6.*(l2.*t5.*t79+l2.*t9.*t78)+m9.*(l2.*t5.*t81+l2.*t9.*t80)+l2.^2.*m2.*(1.0./4.0),m3.*(l2.*l3.*t5.*t7.*(1.0./2.0)+l2.*l3.*t9.*t11.*(1.0./2.0))+m5.*(l2.*t5.*t82+l2.*t9.*t83)+m8.*(l2.*t5.*t84+l2.*t9.*t85)+m6.*(l2.*t5.*t87+l2.*t9.*t86)+m4.*(l2.*t5.*t91+l2.*t9.*t90)+m9.*(l2.*t5.*t89+l2.*t9.*t88)+m7.*(l2.*t5.*t93+l2.*t9.*t92),m4.*(l2.*l4.*t5.*t13.*(1.0./2.0)+l2.*l4.*t9.*t19.*(1.0./2.0))+m6.*(l2.*t5.*t95+l2.*t9.*t96)+m9.*(l2.*t5.*t97+l2.*t9.*t98)+m7.*(l2.*t5.*t100+l2.*t9.*t99)+m5.*(l2.*t5.*t102+l2.*t9.*t101)+m8.*(l2.*t5.*t104+l2.*t9.*t103),m5.*(l2.*l5.*t5.*t15.*(1.0./2.0)+l2.*l5.*t9.*t21.*(1.0./2.0))+m7.*(l2.*t5.*t106+l2.*t9.*t107)+m8.*(l2.*t5.*t109+l2.*t9.*t108)+m6.*(l2.*t5.*t111+l2.*t9.*t110)+m9.*(l2.*t5.*t113+l2.*t9.*t112),m6.*(l2.*l6.*t5.*t17.*(1.0./2.0)+l2.*l6.*t9.*t23.*(1.0./2.0))+m8.*(l2.*t5.*t115+l2.*t9.*t116)+m9.*(l2.*t5.*t118+l2.*t9.*t117)+m7.*(l2.*t5.*t120+l2.*t9.*t119),m7.*(l2.*l7.*t5.*t27.*(1.0./2.0)+l2.*l7.*t9.*t25.*(1.0./2.0))+m9.*(l2.*t5.*t122+l2.*t9.*t123)+m8.*(l2.*t5.*t125+l2.*t9.*t124),m8.*(l2.*l8.*t5.*t31.*(1.0./2.0)+l2.*l8.*t9.*t29.*(1.0./2.0))+m9.*(l2.*t5.*t128+l2.*t9.*t127),l2.*l9.*t5.*t41.*(1.0./2.0)+l2.*l9.*t9.*t45.*(1.0./2.0),I3+m3.*(l3.*t7.*t33.*(1.0./2.0)+l3.*t11.*t35.*(1.0./2.0))+m6.*(l3.*t7.*t37+l3.*t11.*t39)+m9.*(l3.*t7.*t43+l3.*t11.*t47)+m4.*(l3.*t7.*t51+l3.*t11.*t49)+m7.*(l3.*t7.*t55+l3.*t11.*t53)+m5.*(l3.*t7.*t63+l3.*t11.*t61)+m8.*(l3.*t7.*t67+l3.*t11.*t65),I3+m4.*(l3.*t7.*t68+l3.*t11.*t69)+m7.*(l3.*t7.*t70+l3.*t11.*t71)+m5.*(l3.*t7.*t73+l3.*t11.*t72)+m8.*(l3.*t7.*t75+l3.*t11.*t74)+m3.*(l3.*t7.*t77.*(1.0./2.0)+l3.*t11.*t76.*(1.0./2.0))+m6.*(l3.*t7.*t79+l3.*t11.*t78)+m9.*(l3.*t7.*t81+l3.*t11.*t80),I3+m3.*(t7.^2.*t94.*(1.0./4.0)+t11.^2.*t94.*(1.0./4.0))+m5.*(l3.*t7.*t82+l3.*t11.*t83)+m8.*(l3.*t7.*t84+l3.*t11.*t85)+m6.*(l3.*t7.*t87+l3.*t11.*t86)+m4.*(l3.*t7.*t91+l3.*t11.*t90)+m9.*(l3.*t7.*t89+l3.*t11.*t88)+m7.*(l3.*t7.*t93+l3.*t11.*t92),m4.*(l3.*l4.*t7.*t13.*(1.0./2.0)+l3.*l4.*t11.*t19.*(1.0./2.0))+m6.*(l3.*t7.*t95+l3.*t11.*t96)+m9.*(l3.*t7.*t97+l3.*t11.*t98)+m7.*(l3.*t7.*t100+l3.*t11.*t99)+m5.*(l3.*t7.*t102+l3.*t11.*t101)+m8.*(l3.*t7.*t104+l3.*t11.*t103),m5.*(l3.*l5.*t7.*t15.*(1.0./2.0)+l3.*l5.*t11.*t21.*(1.0./2.0))+m7.*(l3.*t7.*t106+l3.*t11.*t107)+m8.*(l3.*t7.*t109+l3.*t11.*t108)+m6.*(l3.*t7.*t111+l3.*t11.*t110)+m9.*(l3.*t7.*t113+l3.*t11.*t112),m6.*(l3.*l6.*t7.*t17.*(1.0./2.0)+l3.*l6.*t11.*t23.*(1.0./2.0))+m8.*(l3.*t7.*t115+l3.*t11.*t116)+m9.*(l3.*t7.*t118+l3.*t11.*t117)+m7.*(l3.*t7.*t120+l3.*t11.*t119),m7.*(l3.*l7.*t7.*t27.*(1.0./2.0)+l3.*l7.*t11.*t25.*(1.0./2.0))+m9.*(l3.*t7.*t122+l3.*t11.*t123)+m8.*(l3.*t7.*t125+l3.*t11.*t124),m8.*(l3.*l8.*t7.*t31.*(1.0./2.0)+l3.*l8.*t11.*t29.*(1.0./2.0))+m9.*(l3.*t7.*t128+l3.*t11.*t127),l3.*l9.*t7.*t41.*(1.0./2.0)+l3.*l9.*t11.*t45.*(1.0./2.0),I4+m6.*(l4.*t13.*t37+l4.*t19.*t39)+m9.*(l4.*t13.*t43+l4.*t19.*t47)+m4.*(l4.*t13.*t51.*(1.0./2.0)+l4.*t19.*t49.*(1.0./2.0))+m7.*(l4.*t13.*t55+l4.*t19.*t53)+m5.*(l4.*t13.*t63+l4.*t19.*t61)+m8.*(l4.*t13.*t67+l4.*t19.*t65),I4+m4.*(l4.*t13.*t68.*(1.0./2.0)+l4.*t19.*t69.*(1.0./2.0))+m7.*(l4.*t13.*t70+l4.*t19.*t71)+m5.*(l4.*t13.*t73+l4.*t19.*t72)+m8.*(l4.*t13.*t75+l4.*t19.*t74)+m6.*(l4.*t13.*t79+l4.*t19.*t78)+m9.*(l4.*t13.*t81+l4.*t19.*t80),I4+m5.*(l4.*t13.*t82+l4.*t19.*t83)+m8.*(l4.*t13.*t84+l4.*t19.*t85)+m6.*(l4.*t13.*t87+l4.*t19.*t86)+m9.*(l4.*t13.*t89+l4.*t19.*t88)+m4.*(l4.*t13.*t91.*(1.0./2.0)+l4.*t19.*t90.*(1.0./2.0))+m7.*(l4.*t13.*t93+l4.*t19.*t92),I4+m4.*(t13.^2.*t105.*(1.0./4.0)+t19.^2.*t105.*(1.0./4.0))+m6.*(l4.*t13.*t95+l4.*t19.*t96)+m9.*(l4.*t13.*t97+l4.*t19.*t98)+m7.*(l4.*t13.*t100+l4.*t19.*t99)+m5.*(l4.*t13.*t102+l4.*t19.*t101)+m8.*(l4.*t13.*t104+l4.*t19.*t103),m5.*(l4.*l5.*t13.*t15.*(1.0./2.0)+l4.*l5.*t19.*t21.*(1.0./2.0))+m7.*(l4.*t13.*t106+l4.*t19.*t107)+m8.*(l4.*t13.*t109+l4.*t19.*t108)+m6.*(l4.*t13.*t111+l4.*t19.*t110)+m9.*(l4.*t13.*t113+l4.*t19.*t112),m6.*(l4.*l6.*t13.*t17.*(1.0./2.0)+l4.*l6.*t19.*t23.*(1.0./2.0))+m8.*(l4.*t13.*t115+l4.*t19.*t116)+m9.*(l4.*t13.*t118+l4.*t19.*t117)+m7.*(l4.*t13.*t120+l4.*t19.*t119),m7.*(l4.*l7.*t13.*t27.*(1.0./2.0)+l4.*l7.*t19.*t25.*(1.0./2.0))+m9.*(l4.*t13.*t122+l4.*t19.*t123)+m8.*(l4.*t13.*t125+l4.*t19.*t124),m8.*(l4.*l8.*t13.*t31.*(1.0./2.0)+l4.*l8.*t19.*t29.*(1.0./2.0))+m9.*(l4.*t13.*t128+l4.*t19.*t127),l4.*l9.*t13.*t41.*(1.0./2.0)+l4.*l9.*t19.*t45.*(1.0./2.0),I5+m6.*(l5.*t15.*t37+l5.*t21.*t39)+m9.*(l5.*t15.*t43+l5.*t21.*t47)+m7.*(l5.*t15.*t55+l5.*t21.*t53)+m5.*(l5.*t15.*t63.*(1.0./2.0)+l5.*t21.*t61.*(1.0./2.0))+m8.*(l5.*t15.*t67+l5.*t21.*t65),I5+m7.*(l5.*t15.*t70+l5.*t21.*t71)+m5.*(l5.*t15.*t73.*(1.0./2.0)+l5.*t21.*t72.*(1.0./2.0))+m8.*(l5.*t15.*t75+l5.*t21.*t74)+m6.*(l5.*t15.*t79+l5.*t21.*t78)+m9.*(l5.*t15.*t81+l5.*t21.*t80),I5+m5.*(l5.*t15.*t82.*(1.0./2.0)+l5.*t21.*t83.*(1.0./2.0))+m8.*(l5.*t15.*t84+l5.*t21.*t85)+m6.*(l5.*t15.*t87+l5.*t21.*t86)+m9.*(l5.*t15.*t89+l5.*t21.*t88)+m7.*(l5.*t15.*t93+l5.*t21.*t92),I5+m6.*(l5.*t15.*t95+l5.*t21.*t96)+m9.*(l5.*t15.*t97+l5.*t21.*t98)+m7.*(l5.*t15.*t100+l5.*t21.*t99)+m5.*(l5.*t15.*t102.*(1.0./2.0)+l5.*t21.*t101.*(1.0./2.0))+m8.*(l5.*t15.*t104+l5.*t21.*t103),I5+m5.*(t15.^2.*t114.*(1.0./4.0)+t21.^2.*t114.*(1.0./4.0))+m7.*(l5.*t15.*t106+l5.*t21.*t107)+m8.*(l5.*t15.*t109+l5.*t21.*t108)+m6.*(l5.*t15.*t111+l5.*t21.*t110)+m9.*(l5.*t15.*t113+l5.*t21.*t112),m6.*(l5.*l6.*t15.*t17.*(1.0./2.0)+l5.*l6.*t21.*t23.*(1.0./2.0))+m8.*(l5.*t15.*t115+l5.*t21.*t116)+m9.*(l5.*t15.*t118+l5.*t21.*t117)+m7.*(l5.*t15.*t120+l5.*t21.*t119),m7.*(l5.*l7.*t15.*t27.*(1.0./2.0)+l5.*l7.*t21.*t25.*(1.0./2.0))+m9.*(l5.*t15.*t122+l5.*t21.*t123)+m8.*(l5.*t15.*t125+l5.*t21.*t124),m8.*(l5.*l8.*t15.*t31.*(1.0./2.0)+l5.*l8.*t21.*t29.*(1.0./2.0))+m9.*(l5.*t15.*t128+l5.*t21.*t127),l5.*l9.*t15.*t41.*(1.0./2.0)+l5.*l9.*t21.*t45.*(1.0./2.0),I6+m6.*(l6.*t17.*t37.*(1.0./2.0)+l6.*t23.*t39.*(1.0./2.0))+m9.*(l6.*t17.*t43+l6.*t23.*t47)+m7.*(l6.*t17.*t55+l6.*t23.*t53)+m8.*(l6.*t17.*t67+l6.*t23.*t65),I6+m7.*(l6.*t17.*t70+l6.*t23.*t71)+m8.*(l6.*t17.*t75+l6.*t23.*t74)+m6.*(l6.*t17.*t79.*(1.0./2.0)+l6.*t23.*t78.*(1.0./2.0))+m9.*(l6.*t17.*t81+l6.*t23.*t80),I6+m8.*(l6.*t17.*t84+l6.*t23.*t85)+m6.*(l6.*t17.*t87.*(1.0./2.0)+l6.*t23.*t86.*(1.0./2.0))+m9.*(l6.*t17.*t89+l6.*t23.*t88)+m7.*(l6.*t17.*t93+l6.*t23.*t92),I6+m6.*(l6.*t17.*t95.*(1.0./2.0)+l6.*t23.*t96.*(1.0./2.0))+m9.*(l6.*t17.*t97+l6.*t23.*t98)+m7.*(l6.*t17.*t100+l6.*t23.*t99)+m8.*(l6.*t17.*t104+l6.*t23.*t103),I6+m7.*(l6.*t17.*t106+l6.*t23.*t107)+m8.*(l6.*t17.*t109+l6.*t23.*t108)+m6.*(l6.*t17.*t111.*(1.0./2.0)+l6.*t23.*t110.*(1.0./2.0))+m9.*(l6.*t17.*t113+l6.*t23.*t112),I6+m6.*(t17.^2.*t121.*(1.0./4.0)+t23.^2.*t121.*(1.0./4.0))+m8.*(l6.*t17.*t115+l6.*t23.*t116)+m9.*(l6.*t17.*t118+l6.*t23.*t117)+m7.*(l6.*t17.*t120+l6.*t23.*t119),m7.*(l6.*l7.*t17.*t27.*(1.0./2.0)+l6.*l7.*t23.*t25.*(1.0./2.0))+m9.*(l6.*t17.*t122+l6.*t23.*t123)+m8.*(l6.*t17.*t125+l6.*t23.*t124),m8.*(l6.*l8.*t17.*t31.*(1.0./2.0)+l6.*l8.*t23.*t29.*(1.0./2.0))+m9.*(l6.*t17.*t128+l6.*t23.*t127),l6.*l9.*t17.*t41.*(1.0./2.0)+l6.*l9.*t23.*t45.*(1.0./2.0),I7+m9.*(l7.*t27.*t43+l7.*t25.*t47)+m7.*(l7.*t25.*t53.*(1.0./2.0)+l7.*t27.*t55.*(1.0./2.0))+m8.*(l7.*t25.*t65+l7.*t27.*t67),I7+m7.*(l7.*t25.*t71.*(1.0./2.0)+l7.*t27.*t70.*(1.0./2.0))+m8.*(l7.*t25.*t74+l7.*t27.*t75)+m9.*(l7.*t25.*t80+l7.*t27.*t81),I7+m8.*(l7.*t25.*t85+l7.*t27.*t84)+m9.*(l7.*t25.*t88+l7.*t27.*t89)+m7.*(l7.*t25.*t92.*(1.0./2.0)+l7.*t27.*t93.*(1.0./2.0)),I7+m9.*(l7.*t25.*t98+l7.*t27.*t97)+m7.*(l7.*t25.*t99.*(1.0./2.0)+l7.*t27.*t100.*(1.0./2.0))+m8.*(l7.*t25.*t103+l7.*t27.*t104),I7+m8.*(l7.*t25.*t108+l7.*t27.*t109)+m7.*(l7.*t25.*t107.*(1.0./2.0)+l7.*t27.*t106.*(1.0./2.0))+m9.*(l7.*t25.*t112+l7.*t27.*t113),I7+m8.*(l7.*t25.*t116+l7.*t27.*t115)+m9.*(l7.*t25.*t117+l7.*t27.*t118)+m7.*(l7.*t25.*t119.*(1.0./2.0)+l7.*t27.*t120.*(1.0./2.0)),I7+m7.*(t25.^2.*t126.*(1.0./4.0)+t27.^2.*t126.*(1.0./4.0))+m9.*(l7.*t25.*t123+l7.*t27.*t122)+m8.*(l7.*t25.*t124+l7.*t27.*t125),m8.*(l7.*l8.*t25.*t29.*(1.0./2.0)+l7.*l8.*t27.*t31.*(1.0./2.0))+m9.*(l7.*t25.*t127+l7.*t27.*t128),l7.*l9.*t27.*t41.*(1.0./2.0)+l7.*l9.*t25.*t45.*(1.0./2.0),I8+m9.*(l8.*t31.*t43+l8.*t29.*t47)+m8.*(l8.*t29.*t65.*(1.0./2.0)+l8.*t31.*t67.*(1.0./2.0)),I8+m8.*(l8.*t29.*t74.*(1.0./2.0)+l8.*t31.*t75.*(1.0./2.0))+m9.*(l8.*t29.*t80+l8.*t31.*t81),I8+m8.*(l8.*t29.*t85.*(1.0./2.0)+l8.*t31.*t84.*(1.0./2.0))+m9.*(l8.*t29.*t88+l8.*t31.*t89),I8+m9.*(l8.*t29.*t98+l8.*t31.*t97)+m8.*(l8.*t29.*t103.*(1.0./2.0)+l8.*t31.*t104.*(1.0./2.0)),I8+m8.*(l8.*t29.*t108.*(1.0./2.0)+l8.*t31.*t109.*(1.0./2.0))+m9.*(l8.*t29.*t112+l8.*t31.*t113),I8+m9.*(l8.*t29.*t117+l8.*t31.*t118)+m8.*(l8.*t29.*t116.*(1.0./2.0)+l8.*t31.*t115.*(1.0./2.0)),I8+m9.*(l8.*t29.*t123+l8.*t31.*t122)+m8.*(l8.*t29.*t124.*(1.0./2.0)+l8.*t31.*t125.*(1.0./2.0)),I8+m8.*(t29.^2.*t129.*(1.0./4.0)+t31.^2.*t129.*(1.0./4.0))+m9.*(l8.*t29.*t127+l8.*t31.*t128),l8.*l9.*t31.*t41.*(1.0./2.0)+l8.*l9.*t29.*t45.*(1.0./2.0),I9+m9.*(l9.*t41.*t43.*(1.0./2.0)+l9.*t45.*t47.*(1.0./2.0)),I9+m9.*(l9.*t41.*t81.*(1.0./2.0)+l9.*t45.*t80.*(1.0./2.0)),I9+m9.*(l9.*t41.*t89.*(1.0./2.0)+l9.*t45.*t88.*(1.0./2.0)),I9+m9.*(l9.*t41.*t97.*(1.0./2.0)+l9.*t45.*t98.*(1.0./2.0)),I9+m9.*(l9.*t41.*t113.*(1.0./2.0)+l9.*t45.*t112.*(1.0./2.0)),I9+m9.*(l9.*t41.*t118.*(1.0./2.0)+l9.*t45.*t117.*(1.0./2.0)),I9+m9.*(l9.*t41.*t122.*(1.0./2.0)+l9.*t45.*t123.*(1.0./2.0)),I9+m9.*(l9.*t41.*t128.*(1.0./2.0)+l9.*t45.*t127.*(1.0./2.0)),I9+t41.^2.*t130.*(1.0./4.0)+t45.^2.*t130.*(1.0./4.0)],[9,9]);