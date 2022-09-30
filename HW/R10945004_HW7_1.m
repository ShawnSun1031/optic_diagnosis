x = 0; y = 0; z = 0; cx = 0; cy = 0; cz = 1; w = 1; us = 480; ua = 37; N = 10000; d = 0.005; g = 0.79;
A = 0; T = 0; R = 0; n0 = 1; n1 = 1.4; n2 = 1.4; us2 = 220; ua2 = 2.2; g2 = 0.79; d2 = 0.2;
deltz = 0.0025; deltr = 0.0025;
wide = 0.3; deep = 0.15;
Arz = zeros(wide/deltr,deep/deltz);
Av = zeros(wide/deltr,deep/deltz);
phirz = zeros(wide/deltr,deep/deltz);
R_T = zeros(5,2);

radius_w = 0.05;

for epoch = 1:5
    Arz = zeros(wide/deltr,deep/deltz);
    Av = zeros(wide/deltr,deep/deltz);
    phirz = zeros(wide/deltr,deep/deltz);
    A = 0; R = 0 ; T = 0;
    for i = 1:N
        x = 0; y = 0; z = 0;
        cx = 0; cy = 0; cz = 1;w=1;
        xx = x; yy = y; zz = z; zz = z;
        cxx = cx; cyy = cy; czz = cz;
        in_tissue = true;
        absorb = false;
        is_reflect = false;
        is_scatter = false;
        is_layer = false;
        
        r = (n0-n1)^2/(n0+n1)^2; %%%% First
        R = R + r*w;  
        w = (1-r)*w;
        
        while ~absorb & in_tissue
            if ~is_reflect & ~is_layer
                if 0<=z & z<d
                    s = -log(rand(1))/(ua+us);
                else
                    s = -log(rand(1))/(ua2+us2);
                end
                cos_theta = 1/(2*g)*(1+g^2-((1-g^2)/(1-g+2*g*rand(1)))^2);
                theta = acos(cos_theta);
                phi = 2*pi*rand(1);
                xx = x;
                x = x + s*cx;
                yy = y;
                y = y + s*cy;
                zz = z;
                z = z + s*cz;
                
            end
            

            if (0<z) & (z<(d+d2)) | is_reflect | is_layer         %%%%% in the tissue
                
                if (0<=zz & zz<d & d<=z & z<(d+d2)) & ~is_layer %%% layer1 to layer2
                    is_layer = true; 
                    s1 = (d-zz)/cz;
                    s2 = (ua+us)/(ua2+us2)*(s-s1);
                    x = xx + cx*(s1+s2);
                    y = yy + cy*(s1+s2);               
                    z = zz + cz*(s1+s2);
                    
                elseif d<= zz & zz<(d+d2) & 0<=z & z<d & ~is_layer  %%% layer2 to layer1
                    is_layer = true; 
                    AA(i) = z;
                    AAA(i) = cz;
                    AAAA(i) = d-zz;
                    s1 = (d-zz)/cz;
                    s2 = (ua2+us2)/(ua+us)*(s-s1);
                    if cz>0
                        cz = -cz;      %%%%%%%%% if layer2 --> layer1 --> out ---> layer1  
%                         R = R - w;    
                    end
                    x = xx + cx*(s1+s2);
                    y = yy + cy*(s1+s2);
                    z = zz + cz*(s1+s2);

               %%%%%%%%%%%%%%%%%%%    
                else
                    is_reflect = false; 
                    is_layer = false;
                    if  0<=z & z<d
                        A = A + w*ua/(ua+us);
                        radius = sqrt(x^2+y^2);
                        ir = ceil(radius/deltr-0.5)+1;
                        iz = ceil(z/deltz-0.5)+1;
                        if 1<=ir & ir<=wide/deltr & 1<=iz & iz<=deep/deltz & is_scatter 
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
                                    cz = cos(theta)*czz/abs(czz);
                                    cxx = cx; cyy = cy; czz =cz;
                                else  
                                    cx = sin(theta)/sqrt(1-czz*czz)*(cxx*czz*cos(phi)-cyy*sin(phi)) + cxx*cos(theta);
                                    cy = sin(theta)/sqrt(1-czz*czz)*(cyy*czz*cos(phi)+cxx*sin(phi)) + cyy*cos(theta);
                                    cz = -sin(theta)*cos(phi)*sqrt(1-czz*czz) + czz*cos(theta); 
                                    cxx = cx; cyy = cy; czz =cz;
                                end
                            else
                                w = 0;
                                absorb = true;
                            end
                        else
                            if abs(cz) > 0.99               
                                cx = sin(theta)*cos(phi);
                                cy = sin(theta)*sin(phi);
                                cz = cos(theta)*czz/abs(czz);
                                cxx = cx; cyy = cy; czz =cz;
                            else
                                cx = sin(theta)/sqrt(1-czz*czz)*(cxx*czz*cos(phi)-cyy*sin(phi)) + cxx*cos(theta);
                                cy = sin(theta)/sqrt(1-czz*czz)*(cyy*czz*cos(phi)+cxx*sin(phi)) + cyy*cos(theta);
                                cz = -sin(theta)*cos(phi)*sqrt(1-czz*czz) + czz*cos(theta); 
                                cxx = cx; cyy = cy; czz =cz; 
                            end
                        end
                    elseif d<=z & z<(d+d2)
                        A = A + w*ua2/(ua2+us2);
                        radius = sqrt(x^2+y^2);
                        ir = ceil(radius/deltr-0.5)+1;
                        iz = ceil(z/deltz-0.5)+1;
                        if 1<=ir & ir<=wide/deltr & 1<=iz & iz<=deep/deltz & is_scatter 
                            Arz(ir,iz) = Arz(ir,iz) + w*ua2/(ua2+us2);
                            Vr = (2*ir+1)*pi*deltz*(deltr^2);
                            Av(ir,iz) = Arz(ir,iz)/Vr/N;
                            phirz(ir,iz) = Av(ir,iz)/ua2;
                        end
                        is_scatter = true;

                        w = w*us2/(ua2+us2);
                        if w < 0.001
                            if rand(1)< 1/20
                                w = 20*w;
                                if abs(cz) > 0.9999
                                    cx = sin(theta)*cos(phi);
                                    cy = sin(theta)*sin(phi);
                                    cz = cos(theta)*czz/abs(czz);
                                    cxx = cx; cyy = cy; czz =cz;
                                else
                                    cx = sin(theta)/sqrt(1-czz*czz)*(cxx*czz*cos(phi)-cyy*sin(phi)) + cxx*cos(theta);
                                    cy = sin(theta)/sqrt(1-czz*czz)*(cyy*czz*cos(phi)+cxx*sin(phi)) + cyy*cos(theta);
                                    cz = -sin(theta)*cos(phi)*sqrt(1-czz*czz) + czz*cos(theta); 
                                    cxx = cx; cyy = cy; czz =cz;
                                end
                            else
                                w = 0;
                                absorb = true;
                            end
                        else
                            if abs(cz) > 0.99               
                                cx = sin(theta)*cos(phi);
                                cy = sin(theta)*sin(phi);
                                cz = cos(theta)*czz/abs(czz);
                                cxx = cx; cyy = cy; czz =cz;
                            else
                                cx = sin(theta)/sqrt(1-czz*czz)*(cxx*czz*cos(phi)-cyy*sin(phi)) + cxx*cos(theta);
                                cy = sin(theta)/sqrt(1-czz*czz)*(cyy*czz*cos(phi)+cxx*sin(phi)) + cyy*cos(theta);
                                cz = -sin(theta)*cos(phi)*sqrt(1-czz*czz) + czz*cos(theta); 
                                cxx = cx; cyy = cy; czz =cz;
                            end
                        end
                        
                        
                    end
                end

                
                
                %%%%%%%%%%%%%%%
                 

            else
                if z>=(d+d2)
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
                    z = (d+d2)-(z-(d+d2));             
                    cz = -cz;
                    czz = cz;
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
                    czz = cz;
                    in_tissue = true;
                    is_reflect = true;
                end

            end

        end
    end
    R_T(epoch,:) = [R/N T/N];
end

bar3(log(phirz))
xlabel('mm');
ylabel('mm');
zlabel('Fluence Rate (1/cm2)');
format short
R_T_std = std(R_T);
R_T_mean = mean(R_T);
fprintf('The R(reflectance) is %2.4f +- %2.4fand T(transmittance) is %2.4f +- %2.4f . \n',R_T_mean(1),R_T_std(1),R_T_mean(2),R_T_std(2))