for epoch = 1:5
    x = 0; y = 0; z = 0; cx = 0; cy = 0; cz = 1; w = 1; us = 90; ua = 10; N = 10000; d = 0.02; g = 0;
    A = 0; T = 0; R = 0; n0 = 1; n1 = 1.5;
    d = 1/(ua+us)*20;
    for i = 1:N
        x = 0; y = 0; z = 0;
        cx = 0; cy = 0; cz = 1;w=1;
        cxx = cx; cyy = cy; czz = cz;
        in_tissue = true;
        absorb = false;
        is_reflect = false;

        r = (n0-n1)^2/(n0+n1)^2;  %%%%% First
        R = R + r*w;
        w = (1-r)*w;

        while ~absorb & in_tissue       
            if ~is_reflect
                s = -log(rand(1))/(ua+us);
                theta = acos(2*rand(1)-1);
                phi = 2*pi*rand(1);
                x = x + s*cx;
                y = y + s*cy;
                z = z + s*cz;
            end
            if (0<z) & (z<d) | is_reflect         %%%%% in the tissue
                is_reflect = false;
                A = A + w*ua/(ua+us);
                w = w*us/(ua+us);
                if w < 0.001
                    if rand(1)< 1/20
                        w = 20*w;
                        if abs(cz) > 0.9999
                            cx = sin(theta)*cos(phi);
                            cy = sin(theta)*sin(phi);
                            cz = cos(theta)*cz/abs(cz);
                            cxx = cx; cyy = cy; czz =cz;
                        else
                            cx = sin(theta)/sqrt(1-czz*czz)*(cxx*czz*cos(phi)-cyy*sin(phi)) + cxx*cos(theta);
                            cy = sin(theta)/sqrt(1-czz*czz)*(cyy*czz*cos(phi)+cxx*sin(phi)) + cyy*cos(theta);
                            cz = -sin(theta)*cos(phi)*sqrt(1-czz*czz) + czz*cos(theta); 
                            cxx = cx; cyy = cy; czz =cz;
                        end
                    else
                        absorb = true;
                    end
                else
                    if abs(cz) > 0.99               
                        cx = sin(theta)*cos(phi);
                        cy = sin(theta)*sin(phi);
                        cz = cos(theta)*cz/abs(cz);
                        cxx = cx; cyy = cy; czz =cz;
                    else
                        cx = sin(theta)/sqrt(1-czz*czz)*(cxx*czz*cos(phi)-cyy*sin(phi)) + cxx*cos(theta);
                        cy = sin(theta)/sqrt(1-czz*czz)*(cyy*czz*cos(phi)+cxx*sin(phi)) + cyy*cos(theta);
                        cz = -sin(theta)*cos(phi)*sqrt(1-czz*czz) + czz*cos(theta); 
                        cxx = cx; cyy = cy; czz =cz;
                    end
                end 

            else
                in_tissue = false;
                if z>=d
                    T=T+w;
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
    RR(epoch) = R/N;
end
format short
R_mean = mean(RR)
