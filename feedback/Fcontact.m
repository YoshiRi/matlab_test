function F = Fcontact(x,v)
% 
K = 10000;
C = 10;

F = max(0,-K*x-v*C)*(x<0);
end