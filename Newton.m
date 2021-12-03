clc
clear all

tic

tol = 1E-5;

count = 1;
x0 = [0.01;-0.1];
vet_x{count} = x0;

while count < 1E6
    
    
    H = Qfunc(vet_x{count});
    S = inv(H);

    g = gfunc(vet_x{count});
    g_transposta = g.';

    d = -S*g;
 
    count = count + 1;
    vet_x{count} = Secaoaurea(@func ,d , vet_x{count-1}); 
    
    %disp(vet_x{count});
    
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
% a = 2;
% b = -3;
% c = -3;
% d = 8;
% 
% Q = [ a b ; c d ];
% end

% --------------FUNÇÃO 2--------------

% function y=func(x)
% a = 10;
% b = 1;
% y = (1 + a - b*x(1) - b*x(2))^2 + (b + x(1) + a*x(2) - b*x(1)*x(2))^2;
% end
% 
% function grad = gfunc(x)
% c = 2*x(1)*(x(2)^2) - 4*x(1)*x(2) + 4*x(1) - 20*x(2)^2 + 20*x(2) - 20;
% d = 2*(x(1)^2)*x(2) - 2*(x(1)^2) - 40*x(1)*x(2) + 20*x(1) + 202*x(2) - 2;
% grad = [c ; d];
% end
% 
% function Q = Qfunc(x)
% a = 10;
% b = 1;
% 
% s = -2*a*b + 2*a*x(2) + 2*(b^2)*(x(2)^2) + 2*(b^2) - 4*b*x(2) + 2;
% ss = 4*a*b*x(2) + 2*a*x(1) + 4*(b^2)*x(1)*x(2) - 4*b*x(1);
% sss = -4*a*b*x(2) + 2*a + 4*(b^2)*x(1)*x(2) - 4*b*x(1);
% ssss = 2*(a^2) - 4*a*b*x(1) + 2*(b^2)*(x(1)^2) + 2*(b^2);
% 
% Q = [ s ss ; sss ssss ];
% end



% --------------FUNÇÃO MOLAS--------------

function y=func(x)
a = (sqrt((30 + x(1))^2 + x(2)^2) - 30)^2;

b = (sqrt((30 - x(1))^2 + x(2)^2) - 30)^2;

y = 450*a + 300*b - 360*x(2);
end

function grad=gfunc(x)

a = (900*(x(1) + 30)*(sqrt((x(1) + 30)^2 + x(2)^2) - 30))/sqrt((x(1) + 30)^2 + x(2)^2) - (600*(30 - x(1))*(sqrt((30 - x(1))^2 + x(2)^2) - 30))/sqrt((30 - x(1))^2 + x(2)^2);

b = 60*(-6 + x(2)*(25 - 300/(sqrt(900-60*x(1)+x(1)^2+x(2)^2)) - 450/(sqrt(900+60*x(1)+x(1)^2+x(2)^2))));

grad = ([a;b]);
end

function Q = Qfunc(x)
a = x(1);
b = x(2);

dfdxdx = (150*(2*a - 60)^2)/((a - 30)^2 + b^2) + (225*(2*a + 60)^2)/((a + 30)^2 + b^2) + (600*(((a - 30)^2 + b^2)^(1/2) - 30))/((a - 30)^2 + b^2)^(1/2) + (900*(((a + 30)^2 + b^2)^(1/2) - 30))/((a + 30)^2 + b^2)^(1/2) - (150*(2*a - 60)^2*(((a - 30)^2 + b^2)^(1/2) - 30))/((a - 30)^2 + b^2)^(3/2) - (225*(2*a + 60)^2*(((a + 30)^2 + b^2)^(1/2) - 30))/((a + 30)^2 + b^2)^(3/2);
dfdxdz = (300*b*(2*a - 60))/((a - 30)^2 + b^2) + (450*b*(2*a + 60))/((a + 30)^2 + b^2) - (300*b*(2*a - 60)*(((a - 30)^2 + b^2)^(1/2) - 30))/((a - 30)^2 + b^2)^(3/2) - (450*b*(2*a + 60)*(((a + 30)^2 + b^2)^(1/2) - 30))/((a + 30)^2 + b^2)^(3/2);
dfdzdx = (300*b*(2*a - 60))/((a - 30)^2 + b^2) + (450*b*(2*a + 60))/((a + 30)^2 + b^2) - (300*b*(2*a - 60)*(((a - 30)^2 + b^2)^(1/2) - 30))/((a - 30)^2 + b^2)^(3/2) - (450*b*(2*a + 60)*(((a + 30)^2 + b^2)^(1/2) - 30))/((a + 30)^2 + b^2)^(3/2);
dfdzdz = (600*b^2)/((a - 30)^2 + b^2) + (900*b^2)/((a + 30)^2 + b^2) + (600*(((a - 30)^2 + b^2)^(1/2) - 30))/((a - 30)^2 + b^2)^(1/2) + (900*(((a + 30)^2 + b^2)^(1/2) - 30))/((a + 30)^2 + b^2)^(1/2) - (600*b^2*(((a - 30)^2 + b^2)^(1/2) - 30))/((a - 30)^2 + b^2)^(3/2) - (900*b^2*(((a + 30)^2 + b^2)^(1/2) - 30))/((a + 30)^2 + b^2)^(3/2);

Q = [dfdxdx dfdxdz ; dfdzdx dfdzdz ];
end





