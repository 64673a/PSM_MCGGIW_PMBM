function [Bern,lik] = detectionBern(Bern,C,model,F)
lik_=zeros(1,F);
F=F-1;
for f=0:F
    [Bern.GGIW,lik_(f+1)] = updateGGIW(Bern.GGIW,C,model,f);
    lik_(f+1) = lik_(f+1) + log(Bern.r) + log(model.Pd) + log(Bern.w_death(end));
    Bern.r = 1;
    
    Bern.t_death = Bern.t_death(end);
    Bern.w_death = 1;
end
[Bern,lik]=IMM_Fusion(Bern,lik_,model);
end
