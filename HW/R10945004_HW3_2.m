x = 0; y = 0; z = 0; cx = 0; cy = 0; cz = 1; w = 1; us = 90; ua = 10; N = 10000; d = 0.02; g = 0.75;
A = 0; T = 0; R = 0;
pie_x = zeros(5,3);
for epoch = 1:5
    for i = 1:N
        x = 0; y = 0; z = 0;
        cx = 0; cy = 0; cz = 1;w=1;
        in_tissue = true;
        absorb = false;
        while ~absorb & in_tissue
            s = -log(rand(1))/(ua+us);
            cos_theta = 1/(2*g)*(1+g^2-((1-g^2)/(1-g+2*g*rand(1)))^2);
            theta = acos(cos_theta);
            phi = 2*pi*rand(1);

            x = x + s*cx;
            y = y + s*cy;
            z = z + s*cz;

            if (0<z) & (z<d)         %%%%% in the tissue  
                A = A + w*ua/(ua+us);
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
                    if abs(cz) > 0.9999               
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
                in_tissue = false;
                if z>=d
                    T=T+w;
                else
                    R=R+w;
                end
            end

        end
    end
    pie_x(epoch,:) = [A T R];
    
end
subplot(3,2,1),pie(pie_x(1,:));
subplot(3,2,2),pie(pie_x(2,:));
subplot(3,2,3),pie(pie_x(3,:));
subplot(3,2,4),pie(pie_x(4,:));
subplot(3,2,5),pie(pie_x(5,:));
legend({'A','T','R'});
