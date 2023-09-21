clear all
close all
A = load('cage_14.dat');
P10 = [A(10,2); A(10,3); A(10,4)]; 
P14 = [A(14,2); A(14,3); A(14,4)];
P5 = [A(5,2); A(5,3); A(5,4)];
P11 = [A(11,2); A(11,3); A(11,4)];
P9 = [A(9,2); A(9,3); A(9,4)]; 
P2 = [A(2,2); A(2,3); A(2,4)];
P12 = [A(12,2); A(12,3); A(12,4)];
P6 = [A(6,2); A(6,3); A(6,4)];

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
            Q5 = k*P5;
            Q2 = k*P2;
            Q9 = k*P9;
            Q10 = k*P10;
            Q11 = k*P11;
            Q14 = k*P14;
            Q6 = k*P6;
            Q12 = k*P12;
            R1 = (Q2+Q5)/2;
            R2 = (Q5+Q9)/2;
            R3 = (Q5+Q10)/2;
            R4 = (Q10+Q11)/2;
            R5 = (Q10+Q14)/2;
            R6 = (Q6+Q11)/2;
            R7 = (Q11+Q12)/2;

            V1 = P2-P5;
            V2 = R1-(P2+P5)/2;
            ang = -theta1;
            [Ro,t]=AxelRot(ang,V1,(P2+P5)/2);
            R1 = Ro*(V2)+(P2+P5)/2;

            V1 = P9-P5;
            V2 = R2-(P9+P5)/2;
            ang = -theta1;
            [Ro,t]=AxelRot(ang,V1,(P9+P5)/2);
            R2 = Ro*(V2)+(P9+P5)/2;

            V1 = P10-P5;
            V2 = R3-(P10+P5)/2;
            ang = -theta1;
            [Ro,t]=AxelRot(ang,V1,(P10+P5)/2);
            R3 = Ro*(V2)+(P10+P5)/2;

            V1 = P11-P10;
            V2 = R4-(P11+P10)/2;
            ang = -theta2;
            [Ro,t]=AxelRot(ang,V1,(P11+P10)/2);
            R4 = Ro*(V2)+(P11+P10)/2;

            V1 = P14-P10;
            V2 = R5-(P14+P10)/2;
            ang = theta2;
            [Ro,t]=AxelRot(ang,V1,(P14+P10)/2);
            R5 = Ro*(V2)+(P14+P10)/2;

            V1 = P12-P11;
            V2 = R7-(P12+P11)/2;
            ang = -theta2;
            [Ro,t]=AxelRot(ang,V1,(P12+P11)/2);
            R7 = Ro*(V2)+(P12+P11)/2;

            x1 = acos(dot(R1-P5,R2-P5)/norm(R1-P5)/norm(R2-P5))*180/pi;
            x2 = acos(dot(R3-P5,R2-P5)/norm(R3-P5)/norm(R2-P5))*180/pi;
            x3 = acos(dot(R3-P5,R1-P5)/norm(R3-P5)/norm(R1-P5))*180/pi;
            x_5 = x1+x2+x3;
            y1 = acos(dot(R3-P10,R4-P10)/norm(R3-P10)/norm(R4-P10))*180/pi;
            y2 = acos(dot(R4-P10,R5-P10)/norm(R4-P10)/norm(R5-P10))*180/pi;
            y3 = acos(dot(R5-P10,R3-P10)/norm(R5-P10)/norm(R3-P10))*180/pi;
            x_10 = y1+y2+y3;
            y1 = acos(dot(R7-P11,R4-P11)/norm(R7-P11)/norm(R4-P11))*180/pi;
            y2 = acos(dot(R4-P11,R6-P11)/norm(R4-P11)/norm(R6-P11))*180/pi;
            y3 = acos(dot(R6-P11,R7-P11)/norm(R6-P11)/norm(R7-P11))*180/pi;
            x_11 = y1+y2+y3;
            if (x_5 >= 278) && (x_5 <= 288) && (x_10 >= 278) && (x_10 <= 288) && (x_11 >= 278) && (x_11 <= 288)
                D1 = acos(dot(R1-P5,R1-P2)/norm(R1-P5)/norm(R1-P2))*180/pi;
                D2 = acos(dot(R4-P10,R4-P11)/norm(R4-P10)/norm(R4-P11))*180/pi;
                aa = (9*D1+12*D2)/21;
                D = [D aa];
                ind = [ind i];
            end
        end
        if isempty(D) == 0
            [M,I] = min(D);
            k = s_fac(ind(I));
            Q5 = k*P5;
            Q2 = k*P2;
            Q9 = k*P9;
            Q10 = k*P10;
            Q11 = k*P11;
            Q14 = k*P14;
            Q6 = k*P6;
            Q12 = k*P12;

            R1 = (Q2+Q5)/2;
            R2 = (Q5+Q9)/2;
            R3 = (Q5+Q10)/2;
            R4 = (Q10+Q11)/2;
            R5 = (Q10+Q14)/2;
            R6 = (Q6+Q11)/2;
            R7 = (Q11+Q12)/2;

            V1 = P2-P5;
            V2 = R1-(P2+P5)/2;
            ang = -theta1;
            [Ro,t]=AxelRot(ang,V1,(P2+P5)/2);
            R1 = Ro*(V2)+(P2+P5)/2;

            V1 = P9-P5;
            V2 = R2-(P9+P5)/2;
            ang = -theta1;
            [Ro,t]=AxelRot(ang,V1,(P9+P5)/2);
            R2 = Ro*(V2)+(P9+P5)/2;

            V1 = P10-P5;
            V2 = R3-(P10+P5)/2;
            ang = -theta1;
            [Ro,t]=AxelRot(ang,V1,(P10+P5)/2);
            R3 = Ro*(V2)+(P10+P5)/2;

            V1 = P11-P10;
            V2 = R4-(P11+P10)/2;
            ang = -theta2;
            [Ro,t]=AxelRot(ang,V1,(P11+P10)/2);
            R4 = Ro*(V2)+(P11+P10)/2;

            V1 = P14-P10;
            V2 = R5-(P14+P10)/2;
            ang = theta2;
            [Ro,t]=AxelRot(ang,V1,(P14+P10)/2);
            R5 = Ro*(V2)+(P14+P10)/2;

            V1 = P12-P11;
            V2 = R7-(P12+P11)/2;
            ang = -theta2;
            [Ro,t]=AxelRot(ang,V1,(P12+P11)/2);
            R7 = Ro*(V2)+(P12+P11)/2;

            D1 = acos(dot(R1-P5,R1-P2)/norm(R1-P5)/norm(R1-P2))*180/pi;
            D2 = acos(dot(R4-P10,R4-P11)/norm(R4-P10)/norm(R4-P11))*180/pi;
            A(jj,kk) = (9*D1+12*D2)/21;
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

% if isempty(D)
%     fprintf('It is not posible');
% else
%     [M,I] = min(D);
%     k = s_fac(ind(I));
%     Q5 = k*P5;
%     Q2 = k*P2;
%     Q9 = k*P9;
%     Q10 = k*P10;
%     Q11 = k*P11;
%     Q14 = k*P14;
%     Q6 = k*P6;
%     Q12 = k*P12;
% 
%     R1 = (Q2+Q5)/2;
%     R2 = (Q5+Q9)/2;
%     R3 = (Q5+Q10)/2;
%     R4 = (Q10+Q11)/2;
%     R5 = (Q10+Q14)/2;
%     R6 = (Q6+Q11)/2;
%     R7 = (Q11+Q12)/2;
% 
%     V1 = P2-P5;
%     V2 = R1-(P2+P5)/2;
%     ang = -theta1;
%     [Ro,t]=AxelRot(ang,V1,(P2+P5)/2);
%     R1 = Ro*(V2)+(P2+P5)/2;
% 
%     V1 = P9-P5;
%     V2 = R2-(P9+P5)/2;
%     ang = -theta1;
%     [Ro,t]=AxelRot(ang,V1,(P9+P5)/2);
%     R2 = Ro*(V2)+(P9+P5)/2;
% 
%     V1 = P10-P5;
%     V2 = R3-(P10+P5)/2;
%     ang = -theta1;
%     [Ro,t]=AxelRot(ang,V1,(P10+P5)/2);
%     R3 = Ro*(V2)+(P10+P5)/2;
% 
%     V1 = P11-P10;
%     V2 = R4-(P11+P10)/2;
%     ang = -theta2;
%     [Ro,t]=AxelRot(ang,V1,(P11+P10)/2);
%     R4 = Ro*(V2)+(P11+P10)/2;
% 
%     V1 = P14-P10;
%     V2 = R5-(P14+P10)/2;
%     ang = theta2;
%     [Ro,t]=AxelRot(ang,V1,(P14+P10)/2);
%     R5 = Ro*(V2)+(P14+P10)/2;
% 
%     V1 = P12-P11;
%     V2 = R7-(P12+P11)/2;
%     ang = -theta2;
%     [Ro,t]=AxelRot(ang,V1,(P12+P11)/2);
%     R7 = Ro*(V2)+(P12+P11)/2;
% 
%     plot3([P5(1) P2(1)], [P5(2) P2(2)], [P5(3) P2(3)], 'k','LineWidth',3)
%     axis equal
%     axis off
%     hold on
%     plot3([P5(1) P10(1)], [P5(2) P10(2)], [P5(3) P10(3)], 'k','LineWidth',3)
%     hold on
%     plot3([P5(1) P9(1)], [P5(2) P9(2)], [P5(3) P9(3)], 'k','LineWidth',3)
%     hold on
%     plot3([P10(1) P5(1)], [P10(2) P5(2)], [P10(3) P5(3)], 'k','LineWidth',3)
%     hold on
%     plot3([P10(1) P11(1)], [P10(2) P11(2)], [P10(3) P11(3)], 'k','LineWidth',3)
%     hold on
%     plot3([P10(1) P14(1)], [P10(2) P14(2)], [P10(3) P14(3)], 'k','LineWidth',3)
%     hold on
%     plot3([P11(1) P6(1)], [P11(2) P6(2)], [P11(3) P6(3)], 'k','LineWidth',3)
%     hold on
%     plot3([P11(1) P12(1)], [P11(2) P12(2)], [P11(3) P12(3)], 'k','LineWidth',3)
%     hold on
% %     plot3([Q5(1) Q2(1)], [Q5(2) Q2(2)], [Q5(3) Q2(3)], ':r','LineWidth',3)
% %     hold on
% %     plot3([Q5(1) Q10(1)], [Q5(2) Q10(2)], [Q5(3) Q10(3)], ':r','LineWidth',3)
% %     hold on
% %     plot3([Q5(1) Q9(1)], [Q5(2) Q9(2)], [Q5(3) Q9(3)], ':r','LineWidth',3)
% %     hold on
% %     plot3([Q10(1) Q5(1)], [Q10(2) Q5(2)], [Q10(3) Q5(3)], ':r','LineWidth',3)
% %     hold on
% %     plot3([Q10(1) Q11(1)], [Q10(2) Q11(2)], [Q10(3) Q11(3)], ':r','LineWidth',3)
% %     hold on
% %     plot3([Q10(1) Q14(1)], [Q10(2) Q14(2)], [Q10(3) Q14(3)], ':r','LineWidth',3)
% %     hold on
% %     plot3([Q11(1) Q6(1)], [Q11(2) Q6(2)], [Q11(3) Q6(3)], ':r','LineWidth',3)
% %     hold on
% %     plot3([Q11(1) Q12(1)], [Q11(2) Q12(2)], [Q11(3) Q12(3)], ':r','LineWidth',3)
% %     hold on
%     plot3([R1(1) P2(1)], [R1(2) P2(2)], [R1(3) P2(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R1(1) P5(1)], [R1(2) P5(2)], [R1(3) P5(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R2(1) P9(1)], [R2(2) P9(2)], [R2(3) P9(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R2(1) P5(1)], [R2(2) P5(2)], [R2(3) P5(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R3(1) P10(1)], [R3(2) P10(2)], [R3(3) P10(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R3(1) P5(1)], [R3(2) P5(2)], [R3(3) P5(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R4(1) P10(1)], [R4(2) P10(2)], [R4(3) P10(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R5(1) P14(1)], [R5(2) P14(2)], [R5(3) P14(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R5(1) P10(1)], [R5(2) P10(2)], [R5(3) P10(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R4(1) P11(1)], [R4(2) P11(2)], [R4(3) P11(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R6(1) P6(1)], [R6(2) P6(2)], [R6(3) P6(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R6(1) P11(1)], [R6(2) P11(2)], [R6(3) P11(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R7(1) P12(1)], [R7(2) P12(2)], [R7(3) P12(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R7(1) P11(1)], [R7(2) P11(2)], [R7(3) P11(3)], 'b','LineWidth',3)
%     hold off
%     x1 = acos(dot(R1-P5,R2-P5)/norm(R1-P5)/norm(R2-P5))*180/pi;
%     x2 = acos(dot(R3-P5,R2-P5)/norm(R3-P5)/norm(R2-P5))*180/pi;
%     x3 = acos(dot(R3-P5,R1-P5)/norm(R3-P5)/norm(R1-P5))*180/pi;
%     x_5 = x1+x2+x3;
%     y1 = acos(dot(R3-P10,R4-P10)/norm(R3-P10)/norm(R4-P10))*180/pi;
%     y2 = acos(dot(R4-P10,R5-P10)/norm(R4-P10)/norm(R5-P10))*180/pi;
%     y3 = acos(dot(R5-P10,R3-P10)/norm(R5-P10)/norm(R3-P10))*180/pi;
%     x_10 = y1+y2+y3;
%     y1 = acos(dot(R7-P11,R4-P11)/norm(R7-P11)/norm(R4-P11))*180/pi;
%     y2 = acos(dot(R4-P11,R6-P11)/norm(R4-P11)/norm(R6-P11))*180/pi;
%     y3 = acos(dot(R6-P11,R7-P11)/norm(R6-P11)/norm(R7-P11))*180/pi;
%     x_11 = y1+y2+y3;
%     x5_10_11 = [x_5 x_10 x_11]
%     D1 = acos(dot(R1-P5,R1-P2)/norm(R1-P5)/norm(R1-P2))*180/pi;
%     D2 = acos(dot(R4-P10,R4-P11)/norm(R4-P10)/norm(R4-P11))*180/pi;
%     deltas = [D1 D2]
%     average_delta = (9*D1+12*D2)/21
% end