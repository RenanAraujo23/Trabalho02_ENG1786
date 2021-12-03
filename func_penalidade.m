function y=func_penalidade(x,penalidade)
    if penalidade == 0
        y = (x(1) - 2)^4 + (x(1) - 2*x(2))^2;
    else
        y = (x(1) - 2)^4 + (x(1) - 2*x(2))^2 + penalidade*((x(1)^2 - x(2))^2);
    end
end