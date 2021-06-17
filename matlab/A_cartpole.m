function A_n = A_cartpole(x, u)
%A_CARTPOLE
%    A_N = A_CARTPOLE(U1,X2,X4)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    13-Dec-2020 22:01:30

t2 = cos(x(2));
t3 = sin(x(2));
t4 = x(4).^2;
t5 = t3.^2;
t6 = t2+t4;
t7 = t5+1.0;
t8 = 1.0./t7;
t9 = t8.^2;
A_n = reshape([0.0,0.0,0.0,0.0,0.0,0.0,-t8.*(t5-t2.*t6)-t2.*t3.*t9.*(u(1)+t3.*t6).*2.0,-t8.*(t2.*2.0-t4.*t5-t3.*u(1)+t2.^2.*t4)+t2.*t3.*t9.*(t3.*2.0+t2.*u(1)+t2.*t3.*t4).*2.0,1.0,0.0,0.0,0.0,0.0,1.0,t3.*t8.*x(4).*2.0,t2.*t3.*t8.*x(4).*-2.0],[4,4]);
