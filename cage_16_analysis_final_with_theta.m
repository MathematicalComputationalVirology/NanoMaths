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
            [Ro,t]=AxelRot(ang,V1,(P10+P15)/2);
            R5 = Ro*(V2)+(P10+P15)/2;

            V1 = P16-P15;
            V2 = R4-(P16+P15)/2;
            ang = theta1;
            [Ro,t]=AxelRot(ang,V1,(P16+P15)/2);
            R4 = Ro*(V2)+(P16+P15)/2;

            V1 = P9-P14;
            V2 = R2-(P9+P14)/2;
            ang = theta2;
            [Ro,t]=AxelRot(ang,V1,(P9+P14)/2);
            R2 = Ro*(V2)+(P9+P14)/2;

            V1 = P13-P14;
            V2 = R1-(P13+P14)/2;
            ang = theta2;
            [Ro,t]=AxelRot(ang,V1,(P13+P14)/2);
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
            [Ro,t]=AxelRot(ang,V1,(P10+P15)/2);
            R5 = Ro*(V2)+(P10+P15)/2;

            V1 = P16-P15;
            V2 = R4-(P16+P15)/2;
            ang = theta1;
            [Ro,t]=AxelRot(ang,V1,(P16+P15)/2);
            R4 = Ro*(V2)+(P16+P15)/2;

            V1 = P9-P14;
            V2 = R2-(P9+P14)/2;
            ang = theta2;
            [Ro,t]=AxelRot(ang,V1,(P9+P14)/2);
            R2 = Ro*(V2)+(P9+P14)/2;

            V1 = P13-P14;
            V2 = R1-(P13+P14)/2;
            ang = theta2;
            [Ro,t]=AxelRot(ang,V1,(P13+P14)/2);
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

% if isempty(D)
%     fprintf('It is not posible');
% else
%     [M,I] = min(D);
%     k = s_fac(ind(I));
%     Q9 = k*P9;
%     Q13 = k*P13;
%     Q14 = k*P14;
%     Q15 = k*P15;
%     Q10 = k*P10;
%     Q16 = k*P16;
%     R1 = (Q13+Q14)/2;
%     R2 = (Q9+Q14)/2;
%     R3 = (Q15+Q14)/2;
%     R4 = (Q16+Q15)/2;
%     R5 = (Q10+Q15)/2;
% 
% 
%     V1 = P10-P15;
%     V2 = R5-(P10+P15)/2;
%     ang = -theta1;
%     [Ro,t]=AxelRot(ang,V1,(P10+P15)/2);
%     R5 = Ro*(V2)+(P10+P15)/2;
% 
%     V1 = P16-P15;
%     V2 = R4-(P16+P15)/2;
%     ang = theta1;
%     [Ro,t]=AxelRot(ang,V1,(P16+P15)/2);
%     R4 = Ro*(V2)+(P16+P15)/2;
% 
%     V1 = P9-P14;
%     V2 = R2-(P9+P14)/2;
%     ang = theta2;
%     [Ro,t]=AxelRot(ang,V1,(P9+P14)/2);
%     R2 = Ro*(V2)+(P9+P14)/2;
% 
%     V1 = P13-P14;
%     V2 = R1-(P13+P14)/2;
%     ang = theta2;
%     [Ro,t]=AxelRot(ang,V1,(P13+P14)/2);
%     R1 = Ro*(V2)+(P13+P14)/2;
% 
%     plot3([P14(1) P9(1)], [P14(2) P9(2)], [P14(3) P9(3)], 'k','LineWidth',3)
%     axis equal
%     axis off
%     hold on
%     plot3([P14(1) P13(1)], [P14(2) P13(2)], [P14(3) P13(3)], 'k','LineWidth',3)
%     hold on
%     plot3([P14(1) P15(1)], [P14(2) P15(2)], [P14(3) P15(3)], 'k','LineWidth',3)
%     hold on
%     plot3([P15(1) P10(1)], [P15(2) P10(2)], [P15(3) P10(3)], 'k','LineWidth',3)
%     hold on
%     plot3([P15(1) P16(1)], [P15(2) P16(2)], [P15(3) P16(3)], 'k','LineWidth',3)
%     hold on
%     plot3([R2(1) P9(1)], [R2(2) P9(2)], [R2(3) P9(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R1(1) P13(1)], [R1(2) P13(2)], [R1(3) P13(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R1(1) P14(1)], [R1(2) P14(2)], [R1(3) P14(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R2(1) P14(1)], [R2(2) P14(2)], [R2(3) P14(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R3(1) P14(1)], [R3(2) P14(2)], [R3(3) P14(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R3(1) P15(1)], [R3(2) P15(2)], [R3(3) P15(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R4(1) P15(1)], [R4(2) P15(2)], [R4(3) P15(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R5(1) P15(1)], [R5(2) P15(2)], [R5(3) P15(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R5(1) P10(1)], [R5(2) P10(2)], [R5(3) P10(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R4(1) P16(1)], [R4(2) P16(2)], [R4(3) P16(3)], 'b','LineWidth',3)
%     hold off
% 
%     x1 = acos(dot(R1-P14,R2-P14)/norm(R1-P14)/norm(R2-P14))*180/pi;
%     x2 = acos(dot(R3-P14,R2-P14)/norm(R3-P14)/norm(R2-P14))*180/pi;
%     x3 = acos(dot(R3-P14,R1-P14)/norm(R3-P14)/norm(R1-P14))*180/pi;
%     x555 = x1+x2+x3;
%     y1 = acos(dot(R3-P15,R4-P15)/norm(R3-P15)/norm(R4-P15))*180/pi;
%     y2 = acos(dot(R4-P15,R5-P15)/norm(R4-P15)/norm(R5-P15))*180/pi;
%     y3 = acos(dot(R5-P15,R3-P15)/norm(R5-P15)/norm(R3-P15))*180/pi;
%     y455 = y1+y2+y3;
%     x555y455 = [x555 y455]
%     D1 = acos(dot(R1-P14,R1-P13)/norm(R1-P14)/norm(R1-P13))*180/pi;
%     D2 = acos(dot(R3-P14,R3-P15)/norm(R3-P14)/norm(R3-P15))*180/pi;
%     D3 = acos(dot(R5-P10,R5-P15)/norm(R5-P10)/norm(R5-P15))*180/pi;
%     deltas = [D1 D2 D3]
%     average_deta = mean(deltas)
% end