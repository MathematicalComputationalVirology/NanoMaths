clear all
close all
A = load('cage_10.dat');
P6 = [A(6,4); A(6,5); A(6,6)]; 
P2 = [A(2,4); A(2,5); A(2,6)];
P3 = [A(3,4); A(3,5); A(3,6)];
P9 = [A(9,4); A(9,5); A(9,6)]; 

n_y = 100;
T_range = linspace(0,10,n_y);
A = zeros(1,n_y);

n_x = 200;
s_fac = linspace(1,2,n_x);
for j = 1:n_y
    theta1 = T_range(j);
    D = [];
    ind = [];
    for i = 1:n_x
        k = s_fac(i);

        Q6 = k*P6;
        Q2 = k*P2;
        Q9 = k*P9;
        Q3 = k*P3;



        R1 = (Q2+Q6)/2;
        R2 = (Q6+Q3)/2;
        R3 = (Q6+Q9)/2;

        V1 = P3-P6;
        V2 = R2-(P3+P6)/2;
        ang = theta1;
        [Ro,t]=AxelRot(ang,V1,(P3+P6)/2);
        R2 = Ro*(V2)+(P3+P6)/2;

        V1 = P9-P6;
        V2 = R3-(P9+P6)/2;
        ang = -theta1;
        [Ro,t]=AxelRot(ang,V1,(P9+P6)/2);
        R3 = Ro*(V2)+(P9+P6)/2;

        x1 = acos(dot(R1-P6,R2-P6)/norm(R1-P6)/norm(R2-P6))*180/pi;
        x2 = acos(dot(R3-P6,R2-P6)/norm(R3-P6)/norm(R2-P6))*180/pi;
        x3 = acos(dot(R3-P6,R1-P6)/norm(R3-P6)/norm(R1-P6))*180/pi;
        x = x1+x2+x3;

        if (x >= 305) && (x <= 312)
            D1 = acos(dot(R1-P6,R1-P2)/norm(R1-P6)/norm(R1-P2))*180/pi;
            D2 = acos(dot(R2-P3,R2-P6)/norm(R2-P3)/norm(R2-P6))*180/pi;
            aa = (D1+2*D2)/3;
            if aa < 90
                D = [D aa];
                ind = [ind i];
            end
        end

    end
    if isempty(D) == 0
        [M,I] = max(D);
        k = s_fac(ind(I));
        Q6 = k*P6;
        Q2 = k*P2;
        Q9 = k*P9;
        Q3 = k*P3;
        R1 = (Q2+Q6)/2;
        R2 = (Q6+Q3)/2;
        R3 = (Q6+Q9)/2;

        V1 = P3-P6;
        V2 = R2-(P3+P6)/2;
        ang = theta1;
        [Ro,t]=AxelRot(ang,V1,(P3+P6)/2);
        R2 = Ro*(V2)+(P3+P6)/2;

        V1 = P9-P6;
        V2 = R3-(P9+P6)/2;
        ang = -theta1;
        [Ro,t]=AxelRot(ang,V1,(P9+P6)/2);
        R3 = Ro*(V2)+(P9+P6)/2;
        D1 = acos(dot(R1-P6,R1-P2)/norm(R1-P6)/norm(R1-P2))*180/pi;
        D2 = acos(dot(R2-P3,R2-P6)/norm(R2-P3)/norm(R2-P6))*180/pi;
        A(j) = (D1+2*D2)/3;
    else
        A(j) = NaN;
    end
end
plot(T_range,A,'k','LineWidth',3)
xlabel('$\theta$','Interpreter','latex','FontSize',20)
ylabel('Average $\delta$','Interpreter','latex','FontSize',20)
% 
% if isempty(D)
%     fprintf('It is not posible');
% else
%     [M,I] = max(D);
%     k = s_fac(ind(I));
%     Q6 = k*P6;
%     Q2 = k*P2;
%     Q9 = k*P9;
%     Q3 = k*P3;
%     R1 = (Q2+Q6)/2;
%     R2 = (Q6+Q3)/2;
%     R3 = (Q6+Q9)/2;
% 
%     V1 = P3-P6;
%     V2 = R2-(P3+P6)/2;
%     ang = theta1;
%     [Ro,t]=AxelRot(ang,V1,(P3+P6)/2);
%     R2 = Ro*(V2)+(P3+P6)/2;
% 
%     V1 = P9-P6;
%     V2 = R3-(P9+P6)/2;
%     ang = -theta1;
%     [Ro,t]=AxelRot(ang,V1,(P9+P6)/2);
%     R3 = Ro*(V2)+(P9+P6)/2;
% 
%     plot3([P6(1) P9(1)], [P6(2) P9(2)], [P6(3) P9(3)], 'k','LineWidth',3)
%     axis equal
%     axis off
%     hold on
%     plot3([P6(1) P3(1)], [P6(2) P3(2)], [P6(3) P3(3)], 'k','LineWidth',3)
%     hold on
%     plot3([P6(1) P2(1)], [P6(2) P2(2)], [P6(3) P2(3)], 'k','LineWidth',3)
%     hold on
% %     plot3([Q6(1) Q9(1)], [Q6(2) Q9(2)], [Q6(3) Q9(3)], ':r','LineWidth',3)
% %     hold on
% %     plot3([Q6(1) Q3(1)], [Q6(2) Q3(2)], [Q6(3) Q3(3)], ':r','LineWidth',3)
% %     hold on
% %     plot3([Q6(1) Q2(1)], [Q6(2) Q2(2)], [Q6(3) Q2(3)], ':r','LineWidth',3)
% %     hold on
%     plot3([R1(1) P2(1)], [R1(2) P2(2)], [R1(3) P2(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R1(1) P6(1)], [R1(2) P6(2)], [R1(3) P6(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R2(1) P3(1)], [R2(2) P3(2)], [R2(3) P3(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R2(1) P6(1)], [R2(2) P6(2)], [R2(3) P6(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R3(1) P9(1)], [R3(2) P9(2)], [R3(3) P9(3)], 'b','LineWidth',3)
%     hold on
%     plot3([R3(1) P6(1)], [R3(2) P6(2)], [R3(3) P6(3)], 'b','LineWidth',3)
%     hold off
%     x1 = acos(dot(R1-P6,R2-P6)/norm(R1-P6)/norm(R2-P6))*180/pi;
%     x2 = acos(dot(R3-P6,R2-P6)/norm(R3-P6)/norm(R2-P6))*180/pi;
%     x3 = acos(dot(R3-P6,R1-P6)/norm(R3-P6)/norm(R1-P6))*180/pi;
%     x = x1+x2+x3
%     D1 = acos(dot(R1-P6,R1-P2)/norm(R1-P6)/norm(R1-P2))*180/pi;
%     D2 = acos(dot(R2-P3,R2-P6)/norm(R2-P3)/norm(R2-P6))*180/pi;
%     deltas = [D1 D2]
%     average_deltas = (D1+2*D2)/3
% end
