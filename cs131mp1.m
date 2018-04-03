# Allan Matthew B. Mariano
# 2015-05804
# CS 131 WFY



################################################
#Part 1: Curve Fitting
################################################

# x values z(m)
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
#Summer Temperature Values T(C)
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



onematrix = [
  1
  1
  1
  1
  1
  1
  1
  1
];


#Gradient (First Derivative of f(x))
function ret = gradient(x, ans)
  ret = (3*ans(4)*x*x)+ (2*ans(3)*x) + ans(2);
endfunction

#A Cholesky Function that does Factorization
  
function ret = cholesky(Tot,p)
  A = Tot'*Tot;
  b = Tot'*p;
  I  = eye(columns(A));
  s = columns(I);



  D = A;
  L = I;

  for j = 1:s-1
     M1 = I;
    for k=j+1:s
       M1(k,j) = -D(k,j)/D(j,j);
   end
    L = M1*L;
    D = M1*D;
  end


  D= diag(diag(D));

  L = inv(L);

  y = L\b;
  x = (D*L')\y;

  ret = x;
endfunction

#convert to vandermonde matrix first
y = x.*x;
z = x.*x.*x;

lmao = horzcat(onematrix,x);
lmao2 = horzcat(lmao,y);
polynomial = horzcat (lmao2,z);

#get a0, a1, a2, a3
printf("1.1 = Matrix which represents the polynomial\n")

polynomial

printf("\n===============================================\n")

printf("1.2 = The coeficients from cholesky LDLT\n")
coefficients = cholesky(polynomial,temp)



printf("a0 = %f\n", coefficients(1))
printf("a1 = %f\n", coefficients(2))
printf("a2 = %f\n", coefficients(3))
printf("a3 = %f\n", coefficients(4))

printf("===============================================\n")



#plotting

j = linspace(0, 28);
for i = 1: 100
  y(i) = polyval(flipud(coefficients), j(i));
end

figure(1)
clf()
hold on;
grid on;
xlabel("Depth (m)");
ylabel("Temperature ('C)");
title("Temperature vs Depth");
plot(x,temp, 'o')
plot(j,y, '-')
xlim([0,30])
ylim([-10,30])
hold off;

# get inflection point
# second derivative of 

printf("1.3 = Get the inflection point\n")

thermocline = (-2*coefficients(3))/(6*coefficients(4));

printf("Thermocline = %f m", thermocline)

printf("\n===============================================\n")

# get gradient (Function above)



a = linspace(0, 30);
for i = 1: 100
  b(i) = gradient(a(i),coefficients);
end

figure(2)
clf()
hold on;
grid on;
xlabel("Depth (m)");
ylabel("Gradient ('C/m)");
title("Gradient vs Depth");
plot(a,b, '-')
ylim([-1,1])
xlim([0,30])
hold off;

printf("1.5 = Get the gradient(First Derivative)\n")
gradient1 = gradient(thermocline, coefficients);

printf("Gradient = %f C/m \n", gradient1)

printf("===============================================\n")

#get heatflux

printf("1.7 = Get the heatflux using the gradient\n")

heatflux = ((-0.01*(gradient1/100))*86400);

printf("Heat flux = %f cal / (cm^2 day)\n", heatflux)







