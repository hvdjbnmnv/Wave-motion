clc
clear
H = 0.129672856039977;%Wave  height
T = 2.6;%Wave period
h = 0.5;%Still water level
D=0.12; %Pile diameter
g=9.81;
%  x1-m(Elliptical parameter),x2-K(Complete elliptic integral of the first
%  kind),x3-L(Wave length),x4-d(Trough depth),x5-E(Complete elliptic integral of the second kind)
equations = @(x) [1+H/x(1)/h*(1-x(1)-x(5)/x(2))+(H/x(1)/h)^2*(-0.5+0.5*x(1)+(0.5-0.25*x(1))*x(5)/x(2))+(H/x(1)/h)^3*(133/200-399/400*x(1)+133/400*x(1)^2+(-233/200+233/200*x(1)-1/25*x(1)^2)*x(5)/x(2)+(0.5-0.25*x(1))*(x(5)/x(2))^2)-x(4)/h;
                 quadl(@(xita) 1 ./ (1 - x(1) * (sin(xita)).^2).^0.5, 0, pi/2) - x(2);
                 quadl(@(xita) (1 - x(1) * (sin(xita)).^2).^0.5, 0, pi/2) - x(5);
                 4*x(2)*((3*H/x(1)/h)^-0.5*(1+H/x(1)/h*(1.25-5/8*x(1)-1.5*(x(5)/x(2))))+(H/x(1)/h)^2*(-15/32+15/32*x(1)-21/128*x(1)^2*(1/8+1/16*x(1))*(x(5)/x(2))+3/8*(x(5)/x(2))^2))-x(3)/h];%最新的文章将-1/16改为1/16
%                  (3*H/x(4)/4/x(1))^0.5*(1+(H/x(4)/x(1))*(0.25-7/8*x(1))+(H/x(4)/x(1))^2*(1/32-11/32*x(1)+111/128*x(1)^2))-2*x(2)*x(4)/x(3)];
% Initial guess
m_guess=1-16*exp(-(3*26/4)^0.5)
K_guess=0.5*log(16/(1-m_guess));
L_guess=(g*h)^0.5*T;
d_guess=0.5;
E_guess=quadl(@(xita) (1 - m_guess^0.5 * (sin(xita)).^2).^0.5, 0, pi/2);
initial_guess = [m_guess, K_guess, L_guess, d_guess,E_guess];

solution = fsolve(equations, initial_guess);

m = solution(1);
K = solution(2);
L = solution(3);
d = solution(4);
E = solution(5);
fprintf('m = %.4f\n', m);
fprintf('K = %.4f\n', K);
fprintf('L = %.4f\n', L);
fprintf('d = %.4f\n', d);
fprintf('E = %.4f\n', E);

[K2,E2] = ellipke(m);
%d=h*(1+H/m/h*(1-m-E2/K2)+(H/m/h)^2*(-0.5+0.5*m+(0.5-0.25*m)*E2/K2)+(H/m/h)^3*(133/200-399/400*m+133/400*m^2+(-233/200+233/200*m-1/25*m^2)*E2/K2+(0.5-0.25*m)*(E2/K2)^2));
%L=h*(4*K2*((3*H/m/h)^-0.5*(1+H/m/h*(1.25-5/8*m-1.5*(E2/K2)))+(H/m/h)^2*(-15/32+15/32*m-21/128*m^2*(1/8-1/16*m)*(E2/K2)+3/8*(E2/K2)^2)));
e=E/K;
c=L/T;
epsilon=H/d;
x=0;
t_range=[0:0.01:10*T]';
alpha=(3*epsilon/4/m)^0.5*(1+(epsilon/m)*(0.25-7/8*m)+(epsilon/m)^2*(1/32-11/32*m+111/128*m^2));
alpha2=2*K/L*d;
delta=4/3*alpha.^2;
for i= 1:size(t_range,1);
    t=t_range(i,1);
    [sn,cn,dn] = ellipj(alpha/d*(x-c*t),m);
    eta3(i,1)=d*(1+epsilon*cn^2+0.75*epsilon^2*cn^2*(cn^2-1)+1/80*epsilon^3*cn^2*(-61/m+111+(61/m-212)*cn^2+101*cn^4));
    z_h(i,:) =[0:0.01:0.37,linspace(0.37, eta3(i,1), 20)]';
    for j =1:size(z_h(i,:),2)
        z=z_h(i,j);
        Y=z;
        U_x1(i,j)=delta*(0.5-m+m*cn.^2);
        U_x2(i,j)=delta.^2*(-19/40+79/40*m-79/40*m^2+cn.^2*(-1.5*m+3*m.^2)-m.^2*cn.^4+...
            +(Y/d).^2*(-0.75*m+0.75*m^2+cn.^2*(1.5*m-3*m.^2)+9/4*m^2*cn.^4));
        U_x3(i,j)=delta^3*(55/112-3471/1120*m+7113/1120*m.^2-2371/560*m.^3+cn.^2*(71/40*m-339/40*m.^2+...
            +339/40*m^3)+cn.^4*(27/10*m^2-27/5*m^3)+6/5*m^3*cn^6+...
            (Y/d)^2*(9/8*m-27/8*m^2+9/4*m^3+cn.^2*(-9/4*m+27/2*m^2-27/2*m^3)+cn.^4*(-75/8*m^2+75/4*m^3)-15/2*m^3*cn.^6)+...
            +(Y/d)^4*(-3/16*m+9/16*m^2-3/8*m^3+cn.^2*(3/8*m-51/16*m^2+51/16*m^3)+cn.^4*(45/16*m^2-45/8*m^3)+45/16*m^3*cn.^6));
        U_x(i,j)=(g*d)^0.5*(U_x1(i,j)+U_x2(i,j)+U_x3(i,j));
    end
end
u_cnoi=  U_x;%Velocity of water particle
