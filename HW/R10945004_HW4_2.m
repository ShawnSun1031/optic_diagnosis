%%
%(a)
n0 = 1 ; n1 = 1.4;
thetaC = asin(n0/n1);
theta1 = 0:0.0001:thetaC;
theta0 = asin(n1*sin(theta1)/n0);
Rs = ((n1*cos(theta1)-n0*cos(theta0))./(n1*cos(theta1)+n0*cos(theta0))).^2;
Rt = ((n1*cos(theta0)-n0*cos(theta1))./(n1*cos(theta0)+n0*cos(theta1))).^2;
r = 1/2*(Rs+Rt);
angle = theta1./pi*180;

fill_angle = thetaC/pi*180:0.01:90;
fill_r = ones(size(fill_angle));

angle = horzcat(angle,fill_angle);
r = horzcat(r,fill_r);

plot(angle,r)
title('Reflection coefficient vs Angle')
xlabel('Angle')
ylabel('Reflection coefficient')


%%
%(b)
%%% tissue to air
n0 = 1.4; nair = 1;
syms theta0
theta1 = asin(n0*sin(theta0)/nair);
thetaC = asin(nair/n0);
r =  1/2*(((n0*cos(theta0)-nair*cos(theta1))/(n0*cos(theta0)+nair*cos(theta1)))^2+((n0*cos(theta1)-nair*cos(theta0))/(n0*cos(theta1)+nair*cos(theta0)))^2);
rd1 = r*sin(2*theta0);
rd2 = sin(2*theta0);
q1 = vpa(int(rd1,0,thetaC));
q2 = vpa(int(rd2,thetaC,0.5*pi));
rd_tissue_to_air = abs(q1 + q2);
fprintf('Reflection coefficient of diffuse light(tissue to air) is %4.4f\n',rd_tissue_to_air)

%%% tissue to water
n0 = 1.4; nair = 1.33;
syms theta0
theta1 = asin(n0*sin(theta0)/nair);
thetaC = asin(nair/n0);
r =  1/2*(((n0*cos(theta0)-nair*cos(theta1))/(n0*cos(theta0)+nair*cos(theta1)))^2+((n0*cos(theta1)-nair*cos(theta0))/(n0*cos(theta1)+nair*cos(theta0)))^2);
rd1 = r*sin(2*theta0);
rd2 = sin(2*theta0);
q1 = vpa(int(rd1,0,thetaC));
q2 = vpa(int(rd2,thetaC,0.5*pi));
rd_tissue_to_water = abs(q1 + q2);
fprintf('Reflection coefficient of diffuse light(tissue to water) is %4.4f\n',rd_tissue_to_water)



