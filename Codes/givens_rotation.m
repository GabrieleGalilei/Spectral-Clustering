function [c, s] = givens_rotation(a, b)
    if b == 0
        c = 1;
        s = 0;
    else
        if abs(b) > abs(a)
            temp = a / b;
            s = 1 / sqrt(1 + temp^2);
            c = s * temp;
        else
            temp = b / a;
            c = 1 / sqrt(1 + temp^2);
            s = c * temp;
        end
    end
end