clearvars

% [x_refract,y_refract]=ginput(8);
M=100;
N=100;
g=10;
refract_prof=@(x,y,M,N) 3*((2^(2/3)+0)-(x/M+y/N).^(2/3));
x_refract =[

    1.5187
    1.7523
   26.5187
   28.6215
   55.4907
   59.2290
   76.2850
   80.4907];
y_refract =[

   70.9112
   60.3972
   77.4533
   63.2009
   81.6589
   70.6776
   87.5000
   79.0888];
hold on
for i=1:4
    x1=x_refract(i*2-1);
    x2=x_refract(i*2);
    y1=y_refract(i*2-1);
    y2=y_refract(i*2);
    plot3([x1,x2],[y1,y2],[1;1],'r','LineWidth',3)
    angle(i)=atan((y2-y1)/(x2-x1));
    h(i)=refract_prof((x1+x2)/2,(y1+y2)/2,M,N);
    c(i)=sqrt(g*h(i));
end
figure
sintheta=sin(-pi/4-angle);
scatter(c,sintheta,'filled')
hold on
ylabel('sin(\theta)')
xlabel('Wave Speed (m/s)')
title('sin(\theta) vs Wave Speed During Refraction')
xlim([0,7.5])
ylim([0,.7])
m=c(:)\sintheta(:);
plot([0,7.5],m*[0,7.5])
RMSE=sqrt(sum((sintheta-m*c).^2)/4)
