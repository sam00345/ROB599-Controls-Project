%%Gives out linearized matrices Al and Bl for nonlinear system at:
%state x(1x6 vector)
%%and input u(1x2 vector)
function[Al,Bl]=linearized_mats(xin,uin)
syms x u y v psvar r delta_f F_x

load('Linearized_expressions.mat','A_lin','B_lin')

Al=double(eval(subs(A_lin,[x u y v psvar r delta_f F_x],[xin uin])))
Bl=double(eval(subs(B_lin,[x u y v psvar r delta_f F_x],[xin uin])))