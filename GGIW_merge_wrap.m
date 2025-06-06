function [w_hat,GGIW_hat] = GGIW_merge_wrap(w,GGIWs)


    [w_hat,GGIW_hat.a,GGIW_hat.b,GGIW_hat.m,GGIW_hat.P,GGIW_hat.v,GGIW_hat.V] = ...
    GGIW_merge(exp(w),[GGIWs.a],[GGIWs.b],[GGIWs.m],...
    reshape([GGIWs.P],[4,4,length(w)]),[GGIWs.v],reshape([GGIWs.V],[2,2,length(w)]));
    GGIW_hat.model=0;
end

