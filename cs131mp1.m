# Allan Matthew B. Mariano
# 2015-05804
# CS 131 WFY



################################################
#Part 1: Curve Fitting
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

#convert to vanermonde matrix

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

function ret = gradient(x, ans)
  ret = (3*ans(4)*x*x)+ (2*ans(3)*x) + ans(2);
endfunction
  
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
lmao3 = horzcat (lmao2,z);

#get a0, a1, a2, a3


ans1 = cholesky(lmao3,temp)


#plotting

j = linspace(0, 28);
for i = 1: 100
  y(i) = polyval(flipud(ans1), j(i));
end

figure(1)
clf()
hold on;
xlabel("Depth (m)");
ylabel("Temperature ('C)");
title("Temperature vs Depth");
plot(x,temp, 'o')
plot(j,y, '-')
xlim([0,30])
ylim([-10,30])
hold off;

# get inflection point
# second derivative of yo ass


inflection = (-2*ans1(3))/(6*ans1(4))


# get gradient (Function above)



a = linspace(0, 30);
for i = 1: 100
  b(i) = gradient(a(i),ans1);
end

figure(2)
clf()
hold on;
xlabel("Depth (m)");
ylabel("Gradient ('C/m)");
title("Gradient vs Depth");
plot(a,b, '-')
ylim([-1,1])
xlim([0,30])
hold off;

gradient1 = gradient(inflection, ans1)

#get heatflux

heatflux = ((-0.01*(gradient1/100))*86400)





