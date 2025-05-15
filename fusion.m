function [Bern,lik]=IMM_Fusion(Bern,idxm,lik,model)
    sumlik=-sum(lik)
    for i=0:length(model.motionmodel.F)-1
        GGIW_.a=Bern{idxm}.GGIW(end-i).a*(-lik(end-i))
    end
    Bern{idxm}.GGIW(end-length(model.motionmodel.F)+1).a
end