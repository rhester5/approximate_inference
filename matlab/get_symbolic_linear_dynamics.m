syms x u real
x = sym('x', [4 1], 'real');
u = sym('u', [1 1], 'real');
f = sym('f', [4 1]);

f(1) = x(3);
f(2) = x(4);
f(3) = (u(1) + sin(x(2)) * (x(4)^2 + cos(x(2))))/(1+sin(x(2))^2);
f(4) = (-u(1)*cos(x(2)) - x(4)^2 * cos(x(2)) * sin(x(2)) - 2*sin(x(2)))/(1+sin(x(2))^2);

A_n = jacobian(f, x);
B_n = jacobian(f, u);
matlabFunction(A_n, 'File', 'A_n_cartpole', 'Vars', {x, u});
matlabFunction(B_n, 'File', 'B_n_cartpole', 'Vars', {x, u});
matlabFunction(f, 'File', 'f_n_cartpole', 'Vars', {x, u});

syms x u k l real
x = sym('x', [12 1], 'real');
u = sym('u', [6 1], 'real');
f = sym('f', [12 1]);

f(1) = x(7);
f(2) = x(8);
f(3) = x(9);
f(4) = x(10);
f(5) = x(11);
f(6) = x(12);
dist21 = sqrt((x(1)-x(3))^2 + (x(2)-x(4))^2);
dist12 = sqrt((x(3)-x(1))^2 + (x(4)-x(2))^2);
dist32 = sqrt((x(3)-x(5))^2 + (x(4)-x(6))^2);
dist23 = sqrt((x(5)-x(3))^2 + (x(6)-x(4))^2);
% angle21 = atan2(x(2)-x(4), x(1)-x(3));
% angle12 = atan2(x(4)-x(2), x(3)-x(1));
% angle32 = atan2(x(4)-x(6), x(3)-x(5));
% angle23 = atan2(x(6)-x(4), x(5)-x(3));
% angle21 = asin((x(2)-x(4))/dist21);
% angle12 = asin((x(4)-x(2))/dist12);
% angle32 = asin((x(4)-x(6))/dist32);
% angle23 = asin((x(6)-x(4))/dist23);
angle21 = atan2(x(2)-x(4), (x(1)-x(3) + 1e-3));
angle12 = atan2(x(4)-x(2), (x(3)-x(1) + 1e-3));
angle32 = atan2(x(4)-x(6), (x(3)-x(5) + 1e-3));
angle23 = atan2(x(6)-x(4), (x(5)-x(3) + 1e-3));
% angle21 = atan((x(2)-x(4))/(x(1)-x(3) + 1e-3));
% angle12 = atan((x(4)-x(2))/(x(3)-x(1) + 1e-3));
% angle32 = atan((x(4)-x(6))/(x(3)-x(5) + 1e-3));
% angle23 = atan((x(6)-x(4))/(x(5)-x(3) + 1e-3));
f(7) = u(1) - k * (dist21 - l) * cos(angle21);
f(8) = u(2) - k * (dist21 - l) * sin(angle21);
f(9) = u(3) - k * (dist12 - l) * cos(angle12) - k * (dist32 - l) * cos(angle32);
f(10) = u(4) - k * (dist12 - l) * sin(angle12) - k * (dist32 - l) * sin(angle32);
f(11) = u(5) - k * (dist23 - l) * cos(angle23);
f(12) = u(6) - k * (dist23 - l) * sin(angle23);

A_n = jacobian(f, x);
B_n = jacobian(f, u);
matlabFunction(A_n, 'File', 'A_n_drone', 'Vars',{x, u, k, l});
matlabFunction(B_n, 'File', 'B_n_drone', 'Vars',{x, u, k, l});
matlabFunction(f, 'File', 'f_n_drone', 'Vars',{x, u, k, l});