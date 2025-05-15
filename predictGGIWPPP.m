function [GGIW] = predictGGIWPPP(GGIW,model)

GGIW.a = GGIW.a/model.eta;
GGIW.b = GGIW.b/model.eta;

F = model.motionmodel.PPPF(GGIW.m);
GGIW.m = F*GGIW.m;
GGIW.P = F*GGIW.P*F' + model.motionmodel.Q;

d = 2;
M = eye(2);
e = exp(-model.Ts/model.tao);
GGIW.v = 2*d + 2 + e*(GGIW.v - 2*d - 2);
GGIW.V = e*M*GGIW.V*M';

end

