function [out] = factorialK(c)
    result = 1;
    n = 1;
    while n < c
        n = n + 1;
        result = result * n;
    end
    out = result;
%end

factorial(3)

end