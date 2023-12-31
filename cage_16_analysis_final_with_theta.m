clear all
close all
A = load('cage_16.dat');
P9 = [A(9,2); A(9,3); A(9,4)]; 
P13 = [A(13,2); A(13,3); A(13,4)]; 
P14 = [A(14,2); A(14,3); A(14,4)]; 
P15 = [A(15,2); A(15,3); A(15,4)]; 
P10 = [A(10,2); A(10,3); A(10,4)]; 
P16 = [A(16,2); A(16,3); A(16,4)]; 

% KWOCA 70 %
S_min = 305;
S_max = 312;

% KWOCA 18 %
% S_min = 278;
% S_max = 288;

tt = cputime;
n_y = 100;
T1_range = linspace(0,20,n_y);
n_z = 100;
T2_range = linspace(0,20,n_z);
A = zeros(n_y,n_z);

n_x = 200;
s_fac = linspace(1,2,n_x);

for jj = 1:n_y
    for kk = 1:n_z
        D = [];
        ind = [];
        theta1 = T1_range(jj);
        theta2 = T2_range(kk);
        for i = 1:n_x
            k = s_fac(i);

            Q9 = k*P9;
            Q13 = k*P13;
            Q14 = k*P14;
            Q15 = k*P15;
            Q10 = k*P10;
            Q16 = k*P16;
            R1 = (Q13+Q14)/2;
            R2 = (Q9+Q14)/2;
            R3 = (Q15+Q14)/2;
            R4 = (Q16+Q15)/2;
            R5 = (Q10+Q15)/2;


            V1 = P10-P15;
            V2 = R5-(P10+P15)/2;
            ang = -theta1;
            u = V1/norm(V1);
            ux = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
            Ro = cosd(ang)*eye(3)+sind(ang)*ux+(1-cosd(ang))*(u*u.');
            R5 = Ro*(V2)+(P10+P15)/2;

            V1 = P16-P15;
            V2 = R4-(P16+P15)/2;
            ang = theta1;
            u = V1/norm(V1);
            ux = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
            Ro = cosd(ang)*eye(3)+sind(ang)*ux+(1-cosd(ang))*(u*u.');
            R4 = Ro*(V2)+(P16+P15)/2;

            V1 = P9-P14;
            V2 = R2-(P9+P14)/2;
            ang = theta2;
            u = V1/norm(V1);
            ux = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
            Ro = cosd(ang)*eye(3)+sind(ang)*ux+(1-cosd(ang))*(u*u.');
            R2 = Ro*(V2)+(P9+P14)/2;

            V1 = P13-P14;
            V2 = R1-(P13+P14)/2;
            ang = theta2;
            u = V1/norm(V1);
            ux = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
            Ro = cosd(ang)*eye(3)+sind(ang)*ux+(1-cosd(ang))*(u*u.');
            R1 = Ro*(V2)+(P13+P14)/2;

            x1 = acos(dot(R1-P14,R2-P14)/norm(R1-P14)/norm(R2-P14))*180/pi;
            x2 = acos(dot(R3-P14,R2-P14)/norm(R3-P14)/norm(R2-P14))*180/pi;
            x3 = acos(dot(R3-P14,R1-P14)/norm(R3-P14)/norm(R1-P14))*180/pi;
            x555 = x1+x2+x3;
            y1 = acos(dot(R3-P15,R4-P15)/norm(R3-P15)/norm(R4-P15))*180/pi;
            y2 = acos(dot(R4-P15,R5-P15)/norm(R4-P15)/norm(R5-P15))*180/pi;
            y3 = acos(dot(R5-P15,R3-P15)/norm(R5-P15)/norm(R3-P15))*180/pi;
            y455 = y1+y2+y3;
            if (x555 >= S_min) && (x555 <= S_max) && (y455 >= S_min) && (y455 <= S_max)
                D1 = acos(dot(R1-P14,R1-P13)/norm(R1-P14)/norm(R1-P13))*180/pi;
                D2 = acos(dot(R3-P14,R3-P15)/norm(R3-P14)/norm(R3-P15))*180/pi;
                D3 = acos(dot(R5-P10,R5-P15)/norm(R5-P10)/norm(R5-P15))*180/pi;
                aa = (D1+D2+D3)/3;
                D = [D aa];
                ind = [ind i];
            end
        end
        if isempty(D) == 0
            [M,I] = min(D);
            k = s_fac(ind(I));
            Q9 = k*P9;
            Q13 = k*P13;
            Q14 = k*P14;
            Q15 = k*P15;
            Q10 = k*P10;
            Q16 = k*P16;
            R1 = (Q13+Q14)/2;
            R2 = (Q9+Q14)/2;
            R3 = (Q15+Q14)/2;
            R4 = (Q16+Q15)/2;
            R5 = (Q10+Q15)/2;


            V1 = P10-P15;
            V2 = R5-(P10+P15)/2;
            ang = -theta1;
            u = V1/norm(V1);
            ux = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
            Ro = cosd(ang)*eye(3)+sind(ang)*ux+(1-cosd(ang))*(u*u.');
            R5 = Ro*(V2)+(P10+P15)/2;

            V1 = P16-P15;
            V2 = R4-(P16+P15)/2;
            ang = theta1;
            u = V1/norm(V1);
            ux = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
            Ro = cosd(ang)*eye(3)+sind(ang)*ux+(1-cosd(ang))*(u*u.');
            R4 = Ro*(V2)+(P16+P15)/2;

            V1 = P9-P14;
            V2 = R2-(P9+P14)/2;
            ang = theta2;
            u = V1/norm(V1);
            ux = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
            Ro = cosd(ang)*eye(3)+sind(ang)*ux+(1-cosd(ang))*(u*u.');
            R2 = Ro*(V2)+(P9+P14)/2;

            V1 = P13-P14;
            V2 = R1-(P13+P14)/2;
            ang = theta2;
            u = V1/norm(V1);
            ux = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
            Ro = cosd(ang)*eye(3)+sind(ang)*ux+(1-cosd(ang))*(u*u.');
            R1 = Ro*(V2)+(P13+P14)/2;

            D1 = acos(dot(R1-P14,R1-P13)/norm(R1-P14)/norm(R1-P13))*180/pi;
            D2 = acos(dot(R3-P14,R3-P15)/norm(R3-P14)/norm(R3-P15))*180/pi;
            D3 = acos(dot(R5-P10,R5-P15)/norm(R5-P10)/norm(R5-P15))*180/pi;
            A(jj,kk) = (D1+D2+D3)/3;
        else
            A(jj,kk) = NaN;
        end
    end
end
surf(T1_range,T2_range,A.','EdgeColor','none')
view(0,90)
colorbar
xlabel('$\theta_1$','Interpreter','latex','FontSize',20)
ylabel('$\theta_2$','Interpreter','latex','FontSize',20)
cputime-tt
