clear

tic
count = 1;
y_resticao = 0;
k = 0;
rp(count)=1;
beta=10;

vet_x{count} = [3;2];

while count < 10E5
 
    % Teste da penalidade
    if func_restricao(vet_x{count}) > 0
        penalidade = rp(count);     
    else        
        penalidade = 0;
    end
            
    disp('Vetor que Minimiza')
    disp(vet_x{count})
    disp('Valor da Função')
    disp(func_x(vet_x{count}))
    disp('Valor da Restrição')
    disp(func_restricao(vet_x{count}))
    disp('------------')
    
    % Minimizando   
    count = count + 1;    
    vet_x{count} = Newton(vet_x{count-1}, penalidade);
    
    % Testar se minimizou a restrição

    if (1/2)*rp(count-1)*(func_restricao(vet_x{count})) <= 0.000001
        disp('----------')
        disp('FIM')
        disp('----------')
        break
    end 
    rp(count)=rp(count-1)*10;   
end

disp('-----------')

disp('Função respeita a restrição?')
if func_restricao(vet_x{count}) <= 0
    disp('Sim')
else
    disp('Não')
end

disp('-----------')

disp('Vetor que Minimiza')
disp(vet_x{count})
disp('Valor da Função')
disp(func_x(vet_x{count}))

disp('-----------')


toc

