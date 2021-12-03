function A=Univariante(x0,penalidade)
    tol = 1E-5;
    count = 1;
    vet_x{count} = x0;

    while count < 1E5

        count = count + 1;
        d = [1;0];
        vet_x{count} = Secaoaurea(d , vet_x{count-1},penalidade);
        %disp(vet_x{count-1});       

        count = count + 1;
        d = [0;1];
        vet_x{count} = Secaoaurea(d , vet_x{count-1},penalidade);
        %disp(vet_x{count-1});

        if abs(func_penalidade(vet_x{count},penalidade) - func_penalidade(vet_x{count-1},penalidade)) < tol
            % disp('Gradiente menor que a tolerância -> FIM')
            A = vet_x{count};
            break
        end

    %     if abs(gfunc(vet_x{count})) < tol
    %         disp('Gradiente menor que a tolerância -> FIM')
    %         break
    %     end

    end
end


