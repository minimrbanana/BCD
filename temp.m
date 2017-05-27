d = 7000; 
EV=(1:d)/sqrt(d);
S = sprandsym(d,7/d,EV);
pause (0.1);
p = symrcm(S);
pause (0.1);
q = colperm(S);
pause (0.1);
r = symamd(S);





