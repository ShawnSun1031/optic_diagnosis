x = 0; y = 0; z = 0; cx = 0; cy = 0; cz = 1; w = 1; us = 414; ua = 6; N = 10000; d = 0.15; g = 0.91;
A = 0; T = 0; R = 0; n0 = 1; n1 = 1.37;
deltz = 0.01; deltr = 0.01;
wide = 0.3; deep = 0.15;
Arz = zeros(wide/deltr,deep/deltz);
Av = zeros(wide/deltr,deep/deltz);
phirz = zeros(wide/deltr,deep/deltz);
R_T = zeros(5,2);
for epoch = 1:5
    Arz = zeros(wide/deltr,deep/deltz);
    Av = zeros(wide/deltr,deep/deltz);
    phirz = zeros(wide/deltr,deep/deltz);
    A = 0; R = 0 ; T = 0;
    for i = 1:N
        x = 0; y = 0; z = 0;
        cx = 0; cy = 0; cz = 1;w=1;
        in_tissue = true;
        absorb = false;
        is_reflect = false;
        is_first = true;
        is_scatter = false;
        while ~absorb & in_tissue
            if is_first
                is_first = false;
                r = (n0-n1)^2/(n0+n1)^2;
                R = R + r*w;  
                w = (1-r)*w;
            end

            if ~is_reflect
                s = -log(rand(1))/(ua+us);
                cos_theta = 1/(2*g)*(1+g^2-((1-g^2)/(1-g+2*g*rand(1)))^2);
                theta = acos(cos_theta);
                phi = 2*pi*rand(1);
                x = x + s*cx;
                y = y + s*cy;
                z = z + s*cz;
            end

            if (0<z) & (z<d) | is_reflect         %%%%% in the tissue
                is_reflect = false;
                A = A + w*ua/(ua+us);
                radius = sqrt(x^2+y^2);
                ir = ceil(radius/deltr-0.5)+1;
                iz = ceil(z/deltz-0.5)+1;
                
                if 1<=ir & ir<=30 & 1<=iz & iz<=15 & is_scatter 
                    Arz(ir,iz) = Arz(ir,iz) + w*ua/(ua+us);
                    Vr = (2*ir+1)*pi*deltz*(deltr^2);
                    Av(ir,iz) = Arz(ir,iz)/Vr/N;
                    phirz(ir,iz) = Av(ir,iz)/ua;
                end
                is_scatter = true;

                w = w*us/(ua+us);
                if w < 0.001
                    if rand(1)< 1/20
                        w = 20*w;
                        if abs(cz) > 0.9999
                            cx = sin(theta)*cos(phi);
                            cy = sin(theta)*sin(phi);
                            cz = cos(theta)*cz/abs(cz);
                        else
                            cx = sin(theta)/sqrt(1-cz*cz)*(cx*cz*cos(phi)-cy*sin(phi)) + cx*cos(theta);
                            cy = sin(theta)/sqrt(1-cz*cz)*(cy*cz*cos(phi)+cx*sin(phi)) + cy*cos(theta);
                            cz = -sin(theta)*cos(phi)*sqrt(1-cz*cz) + cz*cos(theta); 
                        end
                    else
                        w = 0;
                        absorb = true;


                    end
                else
                    if abs(cz) > 0.99               
                        cx = sin(theta)*cos(phi);
                        cy = sin(theta)*sin(phi);
                        cz = cos(theta)*cz/abs(cz);
                    else
                        cx = sin(theta)/sqrt(1-cz*cz)*(cx*cz*cos(phi)-cy*sin(phi)) + cx*cos(theta);
                        cy = sin(theta)/sqrt(1-cz*cz)*(cy*cz*cos(phi)+cx*sin(phi)) + cy*cos(theta);
                        cz = -sin(theta)*cos(phi)*sqrt(1-cz*cz) + cz*cos(theta); 
                    end
                end 

            else
                if z>=d
                    incident_angle = acos(cz);
                    thetaT = asin(n1/n0*sin(incident_angle));
                    thetaC = asin(n0/n1);
                    if incident_angle >= thetaC
                        r=1;
                        t=1-r;
                    else
                        Rs = ((n1*cos(incident_angle)-n0*cos(thetaT))/(n1*cos(incident_angle)+n0*cos(thetaT)))^2;
                        Rt = ((n1*cos(thetaT)-n0*cos(incident_angle))/(n1*cos(thetaT)+n0*cos(incident_angle)))^2;
                        r = 1/2*(Rs+Rt);

    %                     r = 1/2*(sin(incident_angle-thetaT)^2/sin(incident_angle+thetaT)^2+tan(incident_angle-thetaT)^2/tan(incident_angle+thetaT)^2);
                        t = 1-r;
                    end
                    T = T + t*w;
                    w = r*w;
                    z = d-(z-d);
                    cz = -cz;
                    in_tissue = true;
                    is_reflect = true;

                else   
                    incident_angle = acos(-cz);
                    thetaT = asin(n1/n0*sin(incident_angle));
                    thetaC = asin(n0/n1);
                    if incident_angle >= thetaC
                        r=1;
                        t=1-r;
                    else
                        Rs = ((n1*cos(incident_angle)-n0*cos(thetaT))/(n1*cos(incident_angle)+n0*cos(thetaT)))^2;
                        Rt = ((n1*cos(thetaT)-n0*cos(incident_angle))/(n1*cos(thetaT)+n0*cos(incident_angle)))^2;
                        r = 1/2*(Rs+Rt);

    %                     r = 1/2*(sin(incident_angle-thetaT)^2/sin(incident_angle+thetaT)^2+tan(incident_angle-thetaT)^2/tan(incident_angle+thetaT)^2);
                        t = 1-r;
                    end
                    R = R + t*w;
                    w = r*w;
                    z = -z;
                    cz = -cz;
                    in_tissue = true;
                    is_reflect = true;
                end

            end

        end
    end
    R_T(epoch,:) = [R/N T/N];
end

bar3(phirz)
xlabel('mm');
ylabel('mm');
zlabel('# of Photons Absorbed');
format short
R_T_std = std(R_T);
R_T_mean = mean(R_T);
fprintf('The R(reflectance) is %2.4f +- %2.4fand T(transmittance) is %2.4f +- %2.4f . \n',R_T_mean(1),R_T_std(1),R_T_mean(2),R_T_std(2))


