% MK
% BE 3343
% Date: Oct, 2022
% Student ID:
%Single Student Project
% Chapter questions
%%
% Question 1 (chapter 6 question 2)

% Description: in each of the following problems, determine the best
% function y(x) (linear, exponential, or power function) to describe the
% data. Plot the function on the same plot with the data. Label and format
% the plots appropriately.

% a 

% Define x and y 
x1=[25 30 35 40 45]
y1=[5 260 480 745 1100]

% Plot the x and y on rectilinear, loglog and semilogy axes

figure
plot(x1,y1)
title("Rectilenar axes-Linear plot: a")
xlabel("X axis")
ylabel("Y axis")


figure
loglog(x1,y1)
title("Loglog axes-Power plot: a")
xlabel("X axis")
ylabel("Y axis")

figure
semilogy(x1,y1)
title("Semilogy axes-Exponential plot: a")
xlabel("X axis")
ylabel("Y axis")

% Analyze the plots to see which plot gives a straight line or close to a
% straight line

a1=polyfit(x1,y1,1)

% b 

% Define x and y 
x2=[2.5 3 3.5 4 4.5 5 5.5 6 7 8 9 10]
y2=[1500 1220 1050 915 810 745 690 620 520 480 410 390]

% Plot the x and y on rectilinear, loglog and semilogy axes

figure
plot(x2,y2)
title("Rectilenar axes-Linear plot: b")
xlabel("X axis")
ylabel("Y axis")


figure
loglog(x2,y2)
title("Loglog axes-Power plot: b")
xlabel("X axis")
ylabel("Y axis")

figure
semilogy(x2,y2)
title("Semilogy axes-Exponential plot: b")
xlabel("X axis")
ylabel("Y axis")

% Analyze the plots to see which plot gives a straight line or close to a
% straight line

% c

% Define x and y 
x3=[550 600 650 700 750]
y3=[41.2 18.62 8.62 3.92 1.86]

% Plot the x and y on rectilinear, loglog and semilogy axes

figure
plot(x3,y3)
title("Rectilenar axes-Linear plot: c")
xlabel("X axis")
ylabel("Y axis")


figure
loglog(x3,y3)
title("Loglog axes-Power plot: c")
xlabel("X axis")
ylabel("Y axis")

figure
semilogy(x3,y3)
title("Semilogy axes-Exponential plot: c")
xlabel("X axis")
ylabel("Y axis")

% Analyze the plots to see which plot gives a straight line or close to a
% straight line


a=['The data set a fits a rectilinear plot since it forms a straight line and the following linear function can be used y=',num2str(a1(1)),'*x',num2str(a1(2))]
disp(a)
disp("The data set b fits loglog plot because it forms a straight line and the power function needs to be used")
disp("The data set c fits semilogy and loglog plots since it forms straight line on those plots and the exponential function needs to be used")
%%
% Question 2 (chapter 6 question 9)
% Delete all previous information 
clear
clc
close all

% Description: the following data give the drying time T of a certain pint
% as a function of the amount of a certain additive A


%a: Find the first, second, third and fourth degree polynomials that fit
%the data and plot each polynomial with the data. Determine the quality of
%the curve t for each by computing J, S, and r^2

%Define x and y

A=[0 1 2 3 4 5 6 7 8 9] % oz
T=[130 115 110 90 89 89 95 100 110 125] %min

% Plot to observe the original plot of the data

plot(A,T)
title("The original plot of drying time T as a function of additive A amount")
xlabel("X axis")
ylabel("Y axis")

% first-degree polynomial
% use polyfit to obtain the first degree function coefficients
p1=polyfit(A,T,1)
% use polyval to obtain y values when inputting A data to p1 first degree
% polynomial function
x1=linspace(min(A),max(A),100)
y1=polyval(p1,x1)

%plot the original and the first degree polynomial data to analyze the
%plots
figure
plot(A,T,'ko')
hold on
plot(x1,y1,'r')
title("First degree polynomial plotting")
xlabel("X axis")
ylabel("Y axis")
legend("Original plot","First degree polynomial")

% calculate J, S, and r^2 using below formulas
 
for k=1:4;
    coeff1=polyfit(A,T,1);
    J1(k)=sum((polyval(coeff1,A)-T).^2);
    end

mu1=mean(T);
for k=1:4;
    S1(k)=sum((T-mu1).^2);
    r2_1(k)=1-J1(k)/S1(k);    
end
disp("The J value for the first degree polynomial fitting is:")
disp(J1)
disp("The S value for the first degree polynomial fitting is:")
disp(S1)
disp("The r^2 value for the first degree polynomial fitting is:")
disp(r2_1)

% second-degree polynomial
% use polyfit to obtain the second degree function coefficients
p2=polyfit(A,T,2)
% use polyval to obtain y values when inputting A data to p2 second degree
% polynomial function
x2=linspace(min(A),max(A),100)
y2=polyval(p2,x2)

%plot the original and the second degree polynomial data to analyze the
%plots
figure
plot(A,T,'ko')
hold on
plot(x2,y2,'r')
title("Second degree polynomial plotting")
xlabel("X axis")
ylabel("Y axis")
legend("Original plot","Second degree polynomial")

% calculate J, S, and r^2 using below formulas
 
for k=1:4;
    coeff2=polyfit(A,T,2);
    J2(k)=sum((polyval(coeff2,A)-T).^2);
    end

mu2=mean(T);
for k=1:4;
    S2(k)=sum((T-mu2).^2);
    r2_2(k)=1-J2(k)/S2(k);    
end
disp("The J value for the second degree polynomial fitting is:")
disp(J2)
disp("The S value for the second degree polynomial fitting is:")
disp(S2)
disp("The r^2 value for the second degree polynomial fitting is:")
disp(r2_2)

% third-degree polynomial
% use polyfit to obtain the third degree function coefficients
p3=polyfit(A,T,3)
% use polyval to obtain y values when inputting A data to p3 third degree
% polynomial function
x3=linspace(min(A),max(A),100)
y3=polyval(p3,x3)

%plot the original and the third degree polynomial data to analyze the
%plots
figure
plot(A,T,'ko')
hold on
plot(x3,y3,'r')
title("Third degree polynomial plotting")
xlabel("X axis")
ylabel("Y axis")
legend("Original plot","Third degree polynomial")

% calculate J, S, and r^2 using below formulas
 
for k=1:4;
    coeff3=polyfit(A,T,3);
    J3(k)=sum((polyval(coeff3,A)-T).^2);
    end

mu3=mean(T);
for k=1:4;
    S3(k)=sum((T-mu3).^2);
    r2_3(k)=1-J3(k)/S3(k);    
end
disp("The J value for the third degree polynomial fitting is:")
disp(J3)
disp("The S value for the third degree polynomial fitting is:")
disp(S3)
disp("The r^2 value for the third degree polynomial fitting is:")
disp(r2_3)

% fourth-degree polynomial
% use polyfit to obtain the fourth degree function coefficients
p4=polyfit(A,T,4)
% use polyval to obtain y values when inputting A data to p4 fourth degree
% polynomial function
x4=linspace(min(A),max(A),100)
y4=polyval(p4,x4)

%plot the original and the fourth degree polynomial data to analyze the
%plots
figure
plot(A,T,'ko')
hold on
plot(x4,y4,'r')
title("Fourth degree polynomial plotting")
xlabel("X axis")
ylabel("Y axis")
legend("Original plot","Fourth degree polynomial")

% calculate J, S, and r^2 using below formulas
 
for k=1:4;
    coeff4=polyfit(A,T,4);
    J4(k)=sum((polyval(coeff4,A)-T).^2);
    end

mu4=mean(T);
for k=1:4;
    S4(k)=sum((T-mu4).^2);
    r2_4(k)=1-J4(k)/S4(k);    
end
disp("The J value for the fourth degree polynomial fitting is:")
disp(J4)
disp("The S value for the fourth degree polynomial fitting is:")
disp(S4)
disp("The r^2 value for the fourth degree polynomial fitting is:")
disp(r2_4)

disp("After analyzing the r^2 values of first, second, third, and fourth degree polynomials, the r^2 value of fourth polynomial is the closest to 1 thus is the optimal function")

%b: use the polynomial giving the best fit to estimate the amount of
%additive that minimizes the drying time

%plot the original and the fourth degree polynomial data to analyze the
%plots
figure
plot(A,T,'ko')
hold on
plot(x4,y4,'r')
title("Fourth degree polynomial plotting")
xlabel("X axis")
ylabel("Y axis")
legend("Original plot","Fourth degree polynomial")

%obtain the minimum additive A to minimize time T 

[minA minT]=ginput(1)
%%
% Question 3 (chapter 6 question 12)

% Delete all previous information 
clear
clc
close all

% Description:The following representss pressure samplesm in pounds per
% square inch (psi), taken in a fuel line once every second for 10 sec.

% a : fit a first, second, third degree polynomial to these data, Plot the
% curve fits along with the data points.

% define x and y which are t and p, respectively
t=[1 2 3 4 5 6 7 8 9 10]
p=[26.1 27.0 28.2 29.0 29.8 30.6 31.1 31.3 31.0 30.5]


% first-degree polynomial
% use polyfit to obtain the first degree function coefficients
p1=polyfit(t,p,1)
% use polyval to obtain y values when inputting A data to p1 first degree
% polynomial function
x1=linspace(min(t),max(t),100)
y1=polyval(p1,x1)

%plot the original and the first degree polynomial data to analyze the
%plots
figure
plot(t,p,'ko')
hold on
plot(x1,y1,'r')
title("First degree polynomial plotting")
xlabel("X axis")
ylabel("Y axis")
legend("Original plot","First degree polynomial")

% calculate J, S, and r^2 using below formulas
 
for k=1:4;
    coeff1=polyfit(t,p,1);
    J1(k)=sum((polyval(coeff1,t)-p).^2);
    end

mu1=mean(p);
for k=1:4;
    S1(k)=sum((p-mu1).^2);
    r2_1(k)=1-J1(k)/S1(k);    
end
disp("The J value for the first degree polynomial fitting is:")
disp(J1)
disp("The S value for the first degree polynomial fitting is:")
disp(S1)
disp("The r^2 value for the first degree polynomial fitting is:")
disp(r2_1)


% second-degree polynomial
% use polyfit to obtain the second degree function coefficients
p2=polyfit(t,p,2)
% use polyval to obtain y values when inputting A data to p2 second degree
% polynomial function
x2=linspace(min(t),max(t),100)
y2=polyval(p2,x2)

%plot the original and the second degree polynomial data to analyze the
%plots
figure
plot(t,p,'ko')
hold on
plot(x2,y2,'r')
title("Second degree polynomial plotting")
xlabel("X axis")
ylabel("Y axis")
legend("Original plot","Second degree polynomial")

% calculate J, S, and r^2 using below formulas
 
for k=1:4;
    coeff2=polyfit(t,p,2);
    J2(k)=sum((polyval(coeff2,t)-p).^2);
    end

mu2=mean(p);
for k=1:4;
    S2(k)=sum((p-mu2).^2);
    r2_2(k)=1-J2(k)/S2(k);    
end
disp("The J value for the second degree polynomial fitting is:")
disp(J2)
disp("The S value for the second degree polynomial fitting is:")
disp(S2)
disp("The r^2 value for the second degree polynomial fitting is:")
disp(r2_2)



% third-degree polynomial
% use polyfit to obtain the third degree function coefficients
p3=polyfit(t,p,3)
% use polyval to obtain y values when inputting A data to p3 third degree
% polynomial function
x3=linspace(min(t),max(t),100)
y3=polyval(p3,x3)

%plot the original and the first degree polynomial data to analyze the
%plots
figure
plot(t,p,'ko')
hold on
plot(x3,y3,'r')
title("Third degree polynomial plotting")
xlabel("X axis")
ylabel("Y axis")
legend("Original plot","Third degree polynomial")

% calculate J, S, and r^2 using below formulas
 
for k=1:4;
    coeff3=polyfit(t,p,3);
    J3(k)=sum((polyval(coeff3,t)-p).^2);
    end

mu3=mean(p);
for k=1:4;
    S3(k)=sum((p-mu3).^2);
    r2_3(k)=1-J3(k)/S3(k);    
end
disp("The J value for the third degree polynomial fitting is:")
disp(J3)
disp("The S value for the third degree polynomial fitting is:")
disp(S3)
disp("The r^2 value for the third degree polynomial fitting is:")
disp(r2_3)


% b: use the results from part a to predict the pressure at t=11 sec.
% Explain which curve fit gives the most reliable prediction. Consider the
% coefficients of determination and the residuals for each fit in making
% your decision.

disp("The coefficient of determination of the first, second, and third degree polynomial are:")
disp(r2_1)
disp(r2_2)
disp(r2_3)

disp("As can be observed from the above coefficients of determination, third degree polynomials seems to be closest to 1 which is the optimal function")

% To find t=11, use polyval function with p3 third degree polynimal fitted
% function

t_11=polyval(p3,11)

disp("The pressure when t=11 sec is :")
disp(t_11)
%%
% Question 4 (chapter 6 question 14)

% Description: the solubility of salt in water is a funciton of the water
% temperature. Let S represent the solubility of NaCL as grams of salt in
% 100 g of water. Let T be temperature in C. Use the following data to
% obtain a curve for S as a function of T. Usee the fit to estimate S when
% T=25 C

% delete all previous information and entries
clear
clc
close all

% Define x and y which are T and S, respectively
x=[10 20 30 40 50 60 70 80 90]
y=[35 35.6 36.25 36.9 37.5 38.1 38.8 39.4 40]

% the following code will be used to figure out which degree polynomial
% fits the data best

% polyfit is used to find coefficients of first, second, third, fourth, and
% fifth degree polynomials that fits the original data 
p1=polyfit(x,y,1)
p2=polyfit(x,y,2)
p3=polyfit(x,y,3)
p4=polyfit(x,y,4)
p5=polyfit(x,y,5)
x1=linspace(min(x),max(x),100)

% polyval will be used to plot x values on those specific functions

figure
plot(x,y,'ko')
% hold on
% plot(x1,polyval(p1,x1))
% hold on
% plot(x1,polyval(p2,x1))
hold on
plot(x1,polyval(p3,x1))
hold on
plot(x1,polyval(p4,x1))
hold on
plot(x1,polyval(p5,x1))
title("The plots with 1st, 2nd, 3rd, 4th, and 5th degree polynomial fits")
xlabel("X axis")
ylabel("Y axis")
legend("Original data","1st degree","2nd degree","3rd degree","4th degree","5th degree")

% Obtain an equation describing the polynomial fit
p1=polyfit(x,y,1)
m=p1(1)
b=p1(2)

ans=['After observing the plot of 5 polynomial fits, it is obvious that the data fits the first degree polynomial and will have linear function which can be described as follows: S=',num2str(m),'*T+',num2str(b)]
disp(ans)

% find S when T=25 which is finding y when x=25
S=polyval(p1,25)

ansb=['The solubility of NaCl described by S when T is 25 C is ',num2str(S)]
disp(ansb)
%%
% Question 5 (chapter 6 question 16)

% Description: The following function is linear in the parameters a1 and
% a2: y(x)=a1+a2*ln(x). Use the least squares regression with the following
% data to estimate the values of a1 and a2. Use the curve fit to estimate
% the values of y at x=2.5 and at x=11.

% delete all the previous data and entries
clear
clc
close all

% define x and y

x=[1 2 3 4 5 6 7 8 9 10]
y=[10 14 16 18 19 20 21 22 23 23]


% Least squares method from textbook. Important to remember that x=ln(x)

% create symbolic variables a and b
% a is a1 and b is a2
syms a b 
% find the J function and its sum
J=(a+b*log(x)-y).^2
J=sum(J)

% find derivatives with respect to a1 and a2
a1=diff(J,a)
a2=diff(J,b)

% set the integer value to 4 for each value
a1=vpa(a1,4)
a2=vpa(a2,4)

% solve the equations for a and b by setting them equal to 0 
eqns=[a1==0, a2==0];
S=solve(eqns,[a b])

% obtain a and b values which are a and b respectively
a=vpa(S.a,4)
b=vpa(S.b,4)

disp("The a1 was calculated to be: ")
disp(a)

disp("The a2 was calculated to be: ")
disp(b)

% Obtain y values when x=2.5 and x=11
% Important note: x=log(x)
x_25=2.5
y1_25=a+b*log(x_25)

x_11=11
y1_11=a+b*log(x_11)

% Least squares regression method with polyfit which provides coefficients
% for first degree polynomial equation

p1=polyfit(log(x),y,1)

% Obtain values for the x=2.5 and x=11
y2_25=polyval(p1,log(2.5))
y2_11=polyval(p1,log(11))

% Compare the equations and the final y values when x=2.5 and x=11
x1=['Using the least squres method from the textbook, the following equation was obtained: y=',num2str(double(a)),'+',num2str(double(b)),'*ln(x). The y values when x=2.5 and x=11 are ',num2str(double(y1_25)),' and ',num2str(double(y1_11)),',respectively.']
disp(x1)

x2=['Using polyfit which utilizes least square regression method, the following equation was obtained: y=',num2str(p1(2)),'+',num2str(p1(1)),'*ln(x). The y values when x=2.5 and x=11 are ',num2str(y2_25),' and ',num2str(y2_11),',respectively.']
disp(x2)
