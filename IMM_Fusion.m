function [Bern,lik]=IMM_Fusion(Bern,lik_,model)
if isempty(lik_)
    lik=0;
else
    lik=0;
    [lik,likIndex]=max(lik_);
    GGIW_= Bern.GGIW(end-likIndex+1);
%     for i=1:length(lik_)
%         lik=lik+lik_(i);
%     end
%     lik=lik/length(lik_);
%     %likmin=min(lik_)+1;
%     sum_e_lik=0;
%     for i=1:length(lik_)
%         %lik_(i)=lik_(i)-likmin;
%         sum_e_lik=sum_e_lik+exp(lik_(i));
%     end
%     
%     GGIW_.a=0;
%     GGIW_.b=0;
%     GGIW_.m=[0;0;0;0];
%     GGIW_.P=zeros(4,4);
%     GGIW_.v=0;
%     GGIW_.V=0;
%     %GGIW_.model=Bern{idxm}.GGIW(end-idxm+1).model;
%     GGIW_.model=[];
%     for i=0:length(model.motionmodel.F)-1
%         GGIW_.a=GGIW_.a+Bern.GGIW(end-i).a*(exp(lik_(i+1))/sum_e_lik);
%         GGIW_.b=GGIW_.b+Bern.GGIW(end-i).b*(exp(lik_(i+1))/sum_e_lik);
%         GGIW_.m=GGIW_.m+Bern.GGIW(end-i).m.*(exp(lik_(i+1))/sum_e_lik);
%         GGIW_.P=GGIW_.P+Bern.GGIW(end-i).P*(exp(lik_(i+1))/sum_e_lik);
%         GGIW_.v=GGIW_.v+Bern.GGIW(end-i).v*(exp(lik_(i+1))/sum_e_lik);
%         GGIW_.V=GGIW_.V+Bern.GGIW(end-i).V*(exp(lik_(i+1))/sum_e_lik);
%         GGIW_.model=[GGIW_.model,(exp(lik_(i+1))/sum_e_lik)];
%     end
    Bern.GGIW(end-length(model.motionmodel.F)+1)=GGIW_;
    for f=1:length(model.motionmodel.F)-1
        Bern.GGIW(end)=[];
    end
end
end