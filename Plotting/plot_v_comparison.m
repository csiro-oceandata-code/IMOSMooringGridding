clf
col={'k.' 'b.' 'g.' 'm.' 'r.' 'y.' 'c.' 'ko' 'bo' 'go' 'mo' 'ro' 'k.' 'b.' 'g.' 'm.' 'r.' 'y.' 'c.' 'ko' 'bo' 'go' 'mo' 'ro'};
a=1; clear leg;
for iz = 2:length(nameu);
    if ~strcmp(char(nameu(iz-1)),char(nameu(iz)))
        if a==1
            plot(imag(u(:,iz-2)),imag(u(:,iz-1)),char(col(a)))
            hold on
            leg(a)=strcat(nameu(iz-2),' vs ',nameu(iz-1));
            a=a+1;
        end
        plot(imag(u(:,iz-1)),imag(u(:,iz)),char(col(a)))
        leg(a)=strcat(nameu(iz-1),' vs ',nameu(iz));
            
        a=a+1;
        pause
    end
end

grid
legend(leg)
title([moorn ' v'])
