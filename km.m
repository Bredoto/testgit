clc,close all,clear all
% add features fargo
% Number of clusters this is given information
Ncl = 2; 

x1center  = 1;
y1center  = 1;
x2center  = 10;
y2center  = 10;




N = 100;
x = zeros(N,1);
y = zeros(N,1);
lc1 = zeros(N,1);
lc2 = zeros(N,1);



for i=1:N    
    if mod(i,2)==0 
        x(i) = x1center + 2*rand*sin(rand) ;
        y(i) = x1center - 2*rand*sin(rand);
    else
        x(i) = x2center - 2*rand*sin(rand) ;
        y(i) = x2center + 2*rand*sin(rand) ;      
    end
    
end


figure (1)
plot(x,y,'b.')
xlabel('x')
xlabel('y')
grid on


while(true)
xc1 = max(x)*rand;
yc1 = max(y)*rand;
xc2 = max(x)*rand;
yc2 = max(y)*rand;

for i=1:length(x)
    lc1(i) = sqrt( (xc1 - x(i)).^2 + ( yc1 - y(i)).^2  );
    lc2(i) = sqrt( (xc2 - x(i)).^2 + ( yc2 - y(i)).^2  );
end

 if lc1(1)*lc2(1)< 0.5*max(x).^2 break;
 end

end



hold on
plot(xc1,yc1,'ro')
hold on
plot(xc2,yc2,'ro')
xlabel('x')
ylabel('y')


% Calculation the lengths
% for 1 center
lc1 = zeros(N,1);
lc2 = zeros(N,1);

k = 0;
while (true)
k = k+1

for i=1:length(x)
    lc1(i) = sqrt( (xc1 - x(i)).^2 + ( yc1 - y(i)).^2  );
    lc2(i) = sqrt( (xc2 - x(i)).^2 + ( yc2 - y(i)).^2  );
end



[M1,I1] = min(lc1)
[M2,I2] = min(lc2)

xc1_new = x(I1);
xc2_new = x(I2);

yc1_new = y(I1);
yc2_new = y(I2);


if  ((xc1_new - xc1)< 0.01  && (xc2_new - xc2)< 0.01 && (yc1_new - yc1)< 0.1  && (yc2_new - yc2)< 0.1 )  break;
end

xc1 = xc1_new;
xc2 = xc2_new;

yc1 = yc1_new;
yc2 = yc2_new;

end

hold on
plot(xc1,yc1,'go')
hold on
plot(xc2,yc2,'go')
xlabel('x')
ylabel('y')


x1 = zeros(N,1);
y1 = zeros(N,1);
x2 = zeros(N,1);
y2 = zeros(N,1);


% Changing Color For clusters
for i = 1:N
    if (lc1(i)< 2)
        x1(i) = x(i); y1(i) = y(i);
    else
        x2(i) = x(i); y2(i) = y(i);
    end
end

figure(3)
plot(x1,y1,'r.')
grid on
hold on
plot(x2,y2,'g.')
xlabel('x')
ylabel('y')
grid on
hold on
plot(xc1,yc1,'b.')
hold on
plot(xc2,yc2,'b.')

% 2+2 = 10



% add some minor changing
% add minot update with more perfect code

clc,close all,clear all
disp('new features')
% another one
% another 2 codes and upda


% need edd em algorithm



