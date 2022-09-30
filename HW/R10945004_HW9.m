ua_ex = 6.5; us_ex = 412; ua_em = 2; us_em = 323; g = 0.89;
n0 = 1; n1 = 1.4; QY = 1; deltr = 0.01; deltz = 0.01; N = 10000;
Arz = zeros(30,15);
loc = zeros(1,3,N);
A = 0; T = 0; R = 0;
for i = 1:N
    x = 0; y = 0; z = 0; cx = 0; cy = 0; cz = 1; w = 1;
    cxx = cx; cyy = cy; czz = cz;
    absorb = false;
    in_tissue = true;
    is_reflect = false;
    is_first = true;
    r = (n0-n1)^2/(n0+n1)^2;  %%%%% Specular reflection?
    if rand(1) < r
        R = R + w;
        w = 0;
        absorb = true;
    end
    
    while ~absorb & in_tissue
        if ~is_reflect
            s = -log(rand(1))/(ua_ex+us_ex);
            cos_theta = 1/(2*g)*(1+g^2-((1-g^2)/(1-g+2*g*rand(1)))^2);
            theta = acos(cos_theta);
%             theta = acos(2*rand(1)-1);
            phi = 2*pi*rand(1);
            x = x + s*cx;
            y = y + s*cy;
            z = z + s*cz;
        end
        if z>0         %%%%% in the tissue  
            is_reflect = false;
            if rand(1) <= ua_ex/(ua_ex+us_ex)
                A = A + w; 
                absorb = true;
                radius = sqrt(x^2+y^2);
                ir = ceil(radius/deltr-0.5)+1;
                iz = ceil(z/deltz-0.5)+1;
                if 1<=ir & ir <=30 & 1<=iz & iz<=15
                    Arz(ir,iz) = Arz(ir,iz) + w;
                end
                loc(:,:,A) = [x y z];
               
                w = 0;

            else
                 if abs(cz) > 0.99   
                    cx = sin(theta)*cos(phi);
                    cy = sin(theta)*sin(phi);
                    cz = cos(theta)*cz/abs(cz);
                    cxx = cx; cyy = cy; czz = cz;
                else
                    cx = sin(theta)/sqrt(1-czz*czz)*(cxx*czz*cos(phi)-cyy*sin(phi)) + cxx*cos(theta);
                    cy = sin(theta)/sqrt(1-czz*czz)*(cyy*czz*cos(phi)+cxx*sin(phi)) + cyy*cos(theta);
                    cz = -sin(theta)*cos(phi)*sqrt(1-czz*czz) + czz*cos(theta); 
                    cxx = cx; cyy = cy; czz = cz;
                end
            end

        else
            if z <= 0
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
                if rand(1) <= r
                    is_reflect = true;
                    in_tissue = true;
                    z = -z;
                    cz = -cz;
                    czz = cz;
                else
                    in_tissue = false;
                    R = R + w;
                    w = 0;
                end

            end
        end
    end
end
R_em = fluorescence(loc,A);
fprintf('the fraction of excitation photons remitted is %2.4f\n',R/N)
fprintf('the fraction of emission photons remitted is %2.4f\n',R_em/N)
bar3(Arz)
xlabel('mm');
ylabel('mm');
zlabel('# of phtons (1/cm3)');
function R_em = fluorescence(xyz,number)
ua_ex = 6.5; us_ex = 412; ua_em = 2; us_em = 323; g = 0.89;
n0 = 1; n1 = 1.4; QY = 1; deltr = 0.01; deltz = 0.01; 
A_em = 0; R_em = 0; w = 1;
for i = 1:number
    cx = 0; cy = 0; cz = 0; w = 1;
    cxx = cx; cyy = cy; czz = cz;
    x = xyz(1,1,i); y = xyz(1,2,i); z = xyz(1,3,i);
    is_first = true;
    in_tissue = true;
    absorb = false;
    is_reflect = false;
    while ~absorb & in_tissue
        
        if is_first
            is_first = false;
            g = 0.0001;            
        end
        if ~is_reflect
            s = -log(rand(1))/(ua_em+us_em);
            cos_theta = 1/(2*g)*(1+g^2-((1-g^2)/(1-g+2*g*rand(1)))^2);
            g = 0.89;
            theta = acos(cos_theta);
    %             theta = acos(2*rand(1)-1);
            phi = 2*pi*rand(1);
            x = x + s*cx;
            y = y + s*cy;
            z = z + s*cz;
        end
        if (0<z) | is_reflect         %%%%% in the tissue
            is_reflect = false;
            A_em = A_em + w*ua_em/(ua_em+us_em);
            w = w*us_em/(ua_em+us_em);
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
            if z<0
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
                R_em = R_em + t*w;
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
            
end
