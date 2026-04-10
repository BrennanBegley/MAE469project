close all
clear all
clc
%canonical units
mu=1;

a=1;
e=0.01671;
inclination=0.00005*pi/180;
Omega=-11.26064*pi/180;
omega=114.20783*pi/180;
theta=-2.48284*pi/180;

r=a*(1-e^2)/(1+e*cos(theta));

p=a*(1-e^2);

r_pqw=[r*cos(theta) r*sin(theta) 0]'

v_pqw=sqrt(mu/p)*[-sin(theta) (e+cos(theta)) 0]'

R_tilda=[cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(inclination)...
    -cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(inclination)...
    sin(Omega)*sin(inclination);...
    sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(inclination)...
    -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(inclination)...
    -cos(Omega)*sin(inclination);...
    sin(omega)*sin(inclination) cos(omega)*sin(inclination) cos(inclination)]

r_ijk=R_tilda*r_pqw

v_ijk=R_tilda*v_pqw