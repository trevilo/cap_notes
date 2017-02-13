function I = Isquare(tau,Iamp)
t = mod(6*pi*tau,2*pi);
    if t<pi
        I = Iamp;
    else
        I = -Iamp;
    end
end
