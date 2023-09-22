clear all
close all
% coordinates of the cage
A = load('28_cage_1.dat');
c = [];
% number of trimers in the cage
n = 28;
for i = 1:n
    c = [c [A(i,4); A(i,5); A(i,6)]];
end
% scaling factor to get right angles around 3-fold axis 
k = 1.266;
fid = fopen('trimers.dat','w');
arm_len = [];
for i = 1:n
    plot3([c(1,i)],[c(2,i)],[c(3,i)],'ok','MarkerFaceColor','k','MarkerSize',8)
    axis equal
    grid on
    hold on
    for j = 1:n
        if (j ~= i) && (norm(c(:,i)-c(:,j)) < 1.5)%1.39-1.43
            plot3([c(1,i) c(1,j)],[c(2,i) c(2,j)],[c(3,i) c(3,j)],'k','LineWidth',3)
            hold on
            mid = k*(c(:,i)+c(:,j))/2;
            plot3([c(1,i) mid(1)],[c(2,i) mid(2)],[c(3,i) mid(3)],'b','LineWidth',4)
            hold on
            plot3([c(1,j) mid(1)],[c(2,j) mid(2)],[c(3,j) mid(3)],'b','LineWidth',4)
            hold on
            arm1 = norm(c(:,i)-mid);
            arm_len = [arm_len ; arm1];
            arm2 = norm(c(:,j)-mid);
            arm_len = [arm_len ; arm2];
            fprintf(fid,'%g %g %g %g %g %g\n',c(1,i),c(2,i),c(3,i),mid(1),mid(2),mid(3));
        end
    end
end
fclose(fid);

Arm_length_between_min_and_max = (min(arm_len)+max(arm_len))/2

A = load('trimers.dat');
freq = load('freq_angles_I32-10_2.dat');
f = zeros(1,n);
for i = 1:n
    bf = 0;
    c1 = A(3*i-2,1:3);
    c2 = A(3*i-2,4:6);
    c3 = A(3*i-1,4:6);
    c4 = A(3*i,4:6);
    ang1 = acos(dot(c2-c1,c3-c1)/norm(c2-c1)/norm(c3-c1))*180/pi;
    ang2 = acos(dot(c3-c1,c4-c1)/norm(c3-c1)/norm(c4-c1))*180/pi;
    ang3 = acos(dot(c4-c1,c2-c1)/norm(c4-c1)/norm(c2-c1))*180/pi;
    s = ang1+ang2+ang3;
    a1 = floor(ang1);
    a2 = floor(ang2);
    a3 = floor(ang3);
    vec1 = sort([a1,a2,a3]);
    arr = freq(end,1:3)-freq(1,1:3);
    for ii = 1:arr(1)+1
        for jj = 1:arr(2)+1
            for kk = 1:arr(3)+1
                vec2 = sort(freq((ii-1)*33*31+(jj-1)*33+kk,1:3));
                if isequal(vec1,vec2)
                    f(i) = freq((ii-1)*33*31+(jj-1)*33+kk,4);
                    bf = 1;
                    break
                end
            end
            if bf
                break
            end
        end
        if bf
            break
        end
    end
end

mf = min(f)
