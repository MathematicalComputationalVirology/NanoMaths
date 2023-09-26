A = load('26_cage.dat');
c = [];
n = 26;
for i = 1:n
    c = [c [A(i,4); A(i,5); A(i,6)]];
end
l_min = 0.845225;
l_max = 0.975067;

nx = 100;
s_fac = linspace(1.21,1.33,nx);
pos_k = [];
for ii = 1:nx
    k = s_fac(ii);
    arm_len = [];
    for i = 1:n
        for j = 1:n
            if (j ~= i) && (norm(c(:,i)-c(:,j)) < 1.5)
                mid = k*(c(:,i)+c(:,j))/2;
                arm1 = norm(c(:,i)-mid);
                arm_len = [arm_len ; arm1];
                arm2 = norm(c(:,j)-mid);
                arm_len = [arm_len ; arm2];
            end
        end
    end
    if (l_min <= min(arm_len)) && (max(arm_len) <= l_max)
        pos_k = [pos_k k];
    end
end
ll = length(pos_k);
mf = zeros(1,ll);
for rr = 1:ll
    k = pos_k(rr);
    fid = fopen('trimers.dat','w');
    arm_len = [];
    for i = 1:n
        for j = 1:n
            if (j ~= i) && (norm(c(:,i)-c(:,j)) < 1.5)
                mid = k*(c(:,i)+c(:,j))/2;                
                fprintf(fid,'%g %g %g %g %g %g\n',c(1,i),c(2,i),c(3,i),mid(1),mid(2),mid(3));
            end
        end
    end
    fclose(fid);
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
    mf(rr) = min(f);
end
plot(pos_k,mf,'k','LineWidth',3)
xlabel('$k$','Interpreter','latex','FontSize',18)
ylabel('$m_f$','Interpreter','latex','FontSize',18)
