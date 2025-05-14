% Function to calculate metric m for x, given space z and time t. 

function y = calc_metric(x,z,t,m)


if m == 1 % int(int(abs))
    if isnan(x)
        y = NaN;
    else
        y = trapz(z,(trapz(t,abs(x),1)));
    end
    
elseif m == 2 % int(ab(int))
    if isnan(x)
        y = NaN;
    else
        y = trapz(z,abs((trapz(t,x,1))));
    end

end

