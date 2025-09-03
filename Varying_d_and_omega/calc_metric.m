% Function to calculate metric m for x, given space z and time t. 

function y = calc_metric(x,z,t,m)


if m == 1 % int(int(abs)) % net strain or net flux
    if isnan(x)
        y = NaN;
    else
        y = trapz(z,(trapz(t,abs(x),1)));
    end
    
elseif m == 2 % int(ab(int)) % takes absolute value of cumulative strain or flux and integrates over space
    if isnan(x)
        y = NaN;
    else
        y = trapz(z,abs((trapz(t,x,1))));
    end

end

