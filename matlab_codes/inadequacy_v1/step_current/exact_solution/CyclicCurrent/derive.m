syms g I nav

n1 = I/2 + I*g/(1+g) + nav - I/6 +(I/2)*(g/(1+g));

n0 = nav - I/6 +(I/2)*(g/(1+g));

V = ((1+2*g)/(1+g))*n1 - (g/(1+g))*n0 - I*g/((1+g^2));

pretty(V)