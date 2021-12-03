clc
clear all

tic

tol = 1E-5;

count = 1;

% x0 = [2;2];
% x0 = [-1;-3];
% x0 = [10;2];
% x0 = [-2;-3];
% x0 = [0.01;-0.10];

vet_x{count} = x0;

% Inicialização e primeira iteração (k=0)

H = Qfunc(x0);
S = inv(H);

g0 = gfunc(x0);
g0_transposta = g0.';

d0 = -S*g0;
d0_transposta = d0.';

alfa_num = g0_transposta*d0;
alfa_den = d0_transposta*H*d0;
alfa = -alfa_num/alfa_den;

count = count + 1;

x1 = x0 + alfa*d0;

vet_x{count} = x1;
 
while count < 1E5
    
    g1 = gfunc(vet_x{count});
    g1_transposta = g1.';

    delta_x_0 = vet_x{count} - vet_x{count-1};
    delta_x_0_transposta = delta_x_0.';
    delta_g_0 = gfunc(vet_x{count}) - gfunc(vet_x{count-1});

    v1 = gfunc(vet_x{count-1})*(1+alfa*sqrt(abs((d0_transposta*delta_g_0)/(delta_x_0_transposta*gfunc(vet_x{count-1}))))) - gfunc(vet_x{count});
    v1_transposta = v1.';

    w1 = delta_x_0/(delta_x_0_transposta*delta_g_0);
    w1_transposta = w1.';
    
    H = Qfunc(x0);
    S = inv(H);

    d0 = -(eye(2) + w1*v1_transposta)*S*(eye(2) + v1*w1_transposta)*gfunc(vet_x{count});
    d0_transposta = d0.';

    alfa_num = (gfunc(vet_x{count})')*d0;
    alfa_den = d0_transposta*H*d0;
    alfa = -alfa_num/alfa_den;
    
    count = count + 1;

    vet_x{count} =  vet_x{count-1} + alfa*d0;
       
    if abs(gfunc(vet_x{count})) < tol
        disp('Gradiente menor que a tolerância -> FIM')
        break
    end

end

disp('---------------------------');
disp('Ponto de mínimo:');
disp(vet_x{count});
disp('---------------------------');
disp('Valor da função no ponto de mínimo:');
disp(func(vet_x{count}));
disp('---------------------------');
disp('Gradiente da função no ponto de mínimo:');
disp(gfunc(vet_x{count}));
disp('---------------------------');

toc

% --------------FUNÇÃO 1--------------

% function y=func(x)
% y = (x(1).^2) - (3.*x(1).*x(2)) + (4.*x(2).^2) + x(1) - x(2);
% end
% 
% function grad = gfunc(x)
% a = 2*x(1) - 3*x(2) + 1;
% b = -3*x(1) + 8*x(2) - 1;
% grad = [a ; b];
% end
% 
% function Q = Qfunc(x)
% Q = [2 -3 ; -3 8];
% end

% --------------FUNÇÃO 2--------------

function y=func(x)
a = 10;
b = 1;
y = (1 + a - b*x(1) - b*x(2))^2 + (b + x(1) + a*x(2) - b*x(1)*x(2))^2;
end

function grad = gfunc(x)
c = 2*x(1)*(x(2)^2) - 4*x(1)*x(2) + 4*x(1) - 20*x(2)^2 + 20*x(2) - 20;
d = 2*(x(1)^2)*x(2) - 2*(x(1)^2) - 40*x(1)*x(2) + 20*x(1) + 202*x(2) - 2;
grad = [c ; d];
end

function Q = Qfunc(x)
a = 10;
b = 1;

s = -2*a*b + 2*a*x(2) + 2*(b^2)*(x(2)^2) + 2*(b^2) - 4*b*x(2) + 2;
ss = 4*a*b*x(2) + 2*a*x(1) + 4*(b^2)*x(1)*x(2) - 4*b*x(1);
sss = -4*a*b*x(2) + 2*a + 4*(b^2)*x(1)*x(2) - 4*b*x(1);
ssss = 2*(a^2) - 4*a*b*x(1) + 2*(b^2)*(x(1)^2) + 2*(b^2);

Q = [ s ss ; sss ssss ];
end



% --------------FUNÇÃO MOLAS--------------

% function y=func(x)
% a = (sqrt((30 + x(1))^2 + x(2)^2) - 30)^2;
% 
% b = (sqrt((30 - x(1))^2 + x(2)^2) - 30)^2;
% 
% y = 450*a + 300*b - 360*x(2);
% end
% 
% function grad=gfunc(x)
% 
% a = (900*(x(1) + 30)*(sqrt((x(1) + 30)^2 + x(2)^2) - 30))/sqrt((x(1) + 30)^2 + x(2)^2) - (600*(30 - x(1))*(sqrt((30 - x(1))^2 + x(2)^2) - 30))/sqrt((30 - x(1))^2 + x(2)^2);
% 
% b = 60*(-6 + x(2)*(25 - 300/(sqrt(900-60*x(1)+x(1)^2+x(2)^2)) - 450/(sqrt(900+60*x(1)+x(1)^2+x(2)^2))));
% 
% grad = ([a;b]);
% end
% 
% function Q = Qfunc(x)
% a = x(1);
% b = x(2);
% 
% dfdxdx = (150*(2*a - 60)^2)/((a - 30)^2 + b^2) + (225*(2*a + 60)^2)/((a + 30)^2 + b^2) + (600*(((a - 30)^2 + b^2)^(1/2) - 30))/((a - 30)^2 + b^2)^(1/2) + (900*(((a + 30)^2 + b^2)^(1/2) - 30))/((a + 30)^2 + b^2)^(1/2) - (150*(2*a - 60)^2*(((a - 30)^2 + b^2)^(1/2) - 30))/((a - 30)^2 + b^2)^(3/2) - (225*(2*a + 60)^2*(((a + 30)^2 + b^2)^(1/2) - 30))/((a + 30)^2 + b^2)^(3/2);
% dfdxdz = (300*b*(2*a - 60))/((a - 30)^2 + b^2) + (450*b*(2*a + 60))/((a + 30)^2 + b^2) - (300*b*(2*a - 60)*(((a - 30)^2 + b^2)^(1/2) - 30))/((a - 30)^2 + b^2)^(3/2) - (450*b*(2*a + 60)*(((a + 30)^2 + b^2)^(1/2) - 30))/((a + 30)^2 + b^2)^(3/2);
% dfdzdx = (300*b*(2*a - 60))/((a - 30)^2 + b^2) + (450*b*(2*a + 60))/((a + 30)^2 + b^2) - (300*b*(2*a - 60)*(((a - 30)^2 + b^2)^(1/2) - 30))/((a - 30)^2 + b^2)^(3/2) - (450*b*(2*a + 60)*(((a + 30)^2 + b^2)^(1/2) - 30))/((a + 30)^2 + b^2)^(3/2);
% dfdzdz = (600*b^2)/((a - 30)^2 + b^2) + (900*b^2)/((a + 30)^2 + b^2) + (600*(((a - 30)^2 + b^2)^(1/2) - 30))/((a - 30)^2 + b^2)^(1/2) + (900*(((a + 30)^2 + b^2)^(1/2) - 30))/((a + 30)^2 + b^2)^(1/2) - (600*b^2*(((a - 30)^2 + b^2)^(1/2) - 30))/((a - 30)^2 + b^2)^(3/2) - (900*b^2*(((a + 30)^2 + b^2)^(1/2) - 30))/((a + 30)^2 + b^2)^(3/2);
% 
% Q = [dfdxdx dfdxdz ; dfdzdx dfdzdz ];
% end





