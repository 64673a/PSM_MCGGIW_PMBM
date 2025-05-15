function [Bern,lik] = misdetectionBern(Bern,model,F)
lik_=zeros(1,F);
for f=0:F-1
    temp = Bern.w_death(end)*model.Pd*(Bern.GGIW(end-f).b/(Bern.GGIW(end-f).b+1))^Bern.GGIW(end-f).a;
    model.Qd = 1-(1-model.Qd)*Bern.w_death(end);
    
    qD = model.Qd + temp;
    
    w1 = model.Qd/qD;
    w2 = temp/qD;
    
    GGIW1 = Bern.GGIW(end-f);
    GGIW2 = GGIW1;
    GGIW2.b = GGIW2.b + 1;
    
    lik = 1 - Bern.r + Bern.r*qD;
    
    [~,Bern.GGIW(end-f).a,Bern.GGIW(end-f).b] = gammaMerge([w1;w2],[GGIW1.a;GGIW2.a],[GGIW1.b;GGIW2.b]);
    
    Bern.r = Bern.r*qD/lik;
    
    lik_(f+1) = log(lik);
    
    %Updated time of death
    Bern.w_death = [Bern.w_death(1:end-1) Bern.w_death(end)*(qD)]/(1-Bern.w_death(end)*(1-qD));
end
[Bern,lik]=IMM_Fusion(Bern,lik_,model);
end
