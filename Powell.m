function A = Powell(x0, penalidade)
    tol = 1E-5;
    count = 1;
    vet_x{count} = x0;
    
    count = count + 1;
    d = [1;0];
    vet_x{count} = Secaoaurea(d , vet_x{count-1}, penalidade);
    if abs(func_penalidade(vet_x{count},penalidade) - func_penalidade(vet_x{count-1},penalidade)) < tol
        A = vet_x{count};
        return
    end

    count = count + 1;
    d = [0;1];
    vet_x{count} = Secaoaurea(d , vet_x{count-1}, penalidade);

    while count < 1E3
    
        d = vet_x{count} - vet_x{count-2};
        count = count + 1;
        vet_x{count} = Secaoaurea(d, vet_x{count-1}, penalidade);
        
        if abs(func_penalidade(vet_x{count},penalidade) - func_penalidade(vet_x{count-1},penalidade)) < tol
            A = vet_x{count};
            break
        end
    end
end