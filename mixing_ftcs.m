function[tracernew,bottomflux]=mixing_ftcs(tracer,deep,Coeff0,Coeff1,Coeff2)

tracer(end+1,:)=deep;

tracernew = tracer(1:end-1,:).*Coeff0 + ...
    [zeros(1,length(tracer(1,:)));tracer(1:end-2,:)].*Coeff1 + ...
    tracer(2:end,:).*Coeff2;

%bottomflux=(tracer(end,:)-tracer(end-1,:))*BottomCoeff;
