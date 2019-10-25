function [ k ] = kurtsis( x )
%kurtsis º∆À„–≈∫≈«Õ∂»

    if all(x==0)
        k=0;
        return;
    end

    x=x-mean(x);
    e=mean(abs(x).^2);
    if e<eps
        k=0;
        return;
    end
    k=mean(abs(x).^4)/e^2;

    if all(isreal(x))
        k=k-3;
    else
        k=k-2;
    end
end

