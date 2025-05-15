function [Bern,lik] = getLineBern(Line,model)
    center=(Line.start+Line.end)/2;
    GGIW_Line=struct('a',50,'b',5,'m',[center(1);center(2); 0; 0;],'P',diag([10,10,10,10]),'v',50,'V',30*eye(2),'model',0);
    LineW=getWforLine(Line,model);
    wp=log(model.birth.w(1));
    [GGIW_c,lik_c] =updateGGIWforPPP(GGIW_Line,LineW,model);
    w_c = lik_c + wp + log(model.Pd);
    [w_hat,Bern.GGIW] = GGIW_merge_wrap(w_c,GGIW_c);
    M_num=length(LineW);
    if M_num > 1
        Bern.r = 1;
        lik = log(w_hat);
    else
        Bern.r = w_hat/(w_hat+model.lambda_fa);
        lik = log(w_hat+model.lambda_fa);
    end
    
    
end