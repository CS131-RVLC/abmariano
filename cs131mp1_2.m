# Allan Matthew B. Mariano
# 2015-05804
# CS 131 WFY



################################################
#Part 3: Piecewise Cubic Polynomial Interpolation
################################################

x =[
  0
  2.3
  4.9
  9.1
  13.7
  18.3
  22.9
  27.2
];
temp = [
  22.8
  22.8
  22.8
  20.6
  13.9
  11.7
  11.1
  11.1
];

function ret = first(x,a)
  ret = (3*(a(1)*(x^2))) + (2*a(2)*x) + a(3);
endfunction

function ret = second(x,a)
  ret =  6*(a(1))*x + 2*a(2);
endfunction

function ret = third(a)
  ret = 6*a(1);
endfunction

function ret = newton(Polynomial)
  old = 11.612;
  err = old;

  while err > 0.01
   curr = old - (second(old, Polynomial)/third(Polynomial));
   err = abs(curr-old);
    old = curr;
  endwhile

  thermocline = curr;
  
  ret = thermocline
  
  
endfunction


gg =[ 
  x(1)^3 x(1)^2 x(1) 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  x(2)^3 x(2)^2 x(2) 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 x(2)^3 x(2)^2 x(2) 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 x(3)^3 x(3)^2 x(3) 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 x(3)^3 x(3)^2 x(3) 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 x(4)^3 x(4)^2 x(4) 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 x(4)^3 x(4)^2 x(4) 1 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 x(5)^3 x(5)^2 x(5) 1 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(5)^3 x(5)^2 x(5) 1 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(6)^3 x(6)^2 x(6) 1 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(6)^3 x(6)^2 x(6) 1 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(7)^3 x(7)^2 x(7) 1 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(7)^3 x(7)^2 x(7) 1
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x(8)^3 x(8)^2 x(8) 1
  (3*(x(2)^2)) (2*(x(2))) 1 0 -(3*(x(2)^2)) -(2*(x(2))) -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 (3*(x(3)^2)) (2*(x(3))) 1 0 -(3*(x(3)^2)) -(2*(x(3))) -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 (3*(x(4)^2)) (2*(x(4))) 1 0 -(3*(x(4)^2)) -(2*(x(4))) -1 0 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 (3*(x(5)^2)) (2*(x(5))) 1 0 -(3*(x(5)^2)) -(2*(x(5))) -1 0 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 (3*(x(6)^2)) (2*(x(6))) 1 0 -(3*(x(6)^2)) -(2*(x(6))) -1 0 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 (3*(x(7)^2)) (2*(x(7))) 1 0 -(3*(x(7)^2)) -(2*(x(7))) -1 0
  6*x(2) 2 0 0 -6*x(2) -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 6*x(3) 2 0 0 -6*x(3) -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 6*x(4) 2 0 0 -6*x(4) -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 6*x(5) 2 0 0 -6*x(5) -2 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6*x(6) 2 0 0 -6*x(6) -2 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6*x(7) 2 0 0 -6*x(7) -2 0 0
  3*(x(1)**2) 2*x(1) 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3*(x(8)**2) 2*x(8) 1 0 
];

b = [
  22.8
  22.8
  22.8
  22.8
  22.8
  20.6
  20.6
  13.9
  13.9
  11.7
  11.7
  11.1
  11.1
  11.1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
];
  
x = gg \ b




FirstPolynomial = [x(1) x(2) x(3) x(4)];
SecondPolynomial = [x(5) x(6) x(7) x(8)];
ThirdPolynomial = [x(9) x(10) x(11) x(12)];
FourthPolynomial = [x(13) x(14) x(15) x(16)];
FifthPolynomial = [x(17) x(18) x(19) x(20)];
SixthPolynomial = [x(21) x(22) x(23) x(24)];
SeventhPolynomial = [x(25) x(26) x(27) x(28)];

t1 = newton(FirstPolynomial);
t2 = newton(SecondPolynomial);
t3 = newton(ThirdPolynomial);
t4 = newton(FourthPolynomial);
t5 = newton(FifthPolynomial);
t6 = newton(SixthPolynomial);
t7 = newton(SeventhPolynomial);


gradient = first(t1, FirstPolynomial)
gradient = first(t2, SecondPolynomial)
gradient = first(t3, ThirdPolynomial)
gradient = first(t4, FourthPolynomial)
gradient = first(t5, FifthPolynomial)
gradient = first(t6, SixthPolynomial)
gradient = first(t7, SeventhPolynomial)

gradient = first(t4, FourthPolynomial)

heatflux = ((-0.01*(gradient/100))*86400)

#plotting Temperature vs Depth

j = linspace(0, 28);
for i = 1: columns(FourthPolynomial)
  y(i) = polyval(flipud(FourthPolynomial), j(i));
end

figure(1)
clf()
hold on;
xlabel("Depth (m)");
ylabel("Temperature ('C)");
title("Temperature vs Depth");
#plot(x,temp, 'o')
plot(j,y, '-')
xlim([0,30])
ylim([-10,30])
hold off;


#plotting the gradient

a = linspace(0, 30);
for i = 1: 100
  b(i) = first(a(i),FourthPolynomial);
end

figure(2)
clf()
xlabel("Depth (m)");
ylabel("Gradient ('C/m)");
title("Gradient vs Depth");
hold on;
plot(a,b, '-')
ylim([-5,5])
xlim([0,30])
hold off;

