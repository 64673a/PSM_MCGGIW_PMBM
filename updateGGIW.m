function [GGIW,lik] = updateGGIW(GGIW,W,model,f)

d = 2;

card_W = size(W,2);

GGIW_.a = GGIW(end-f).a + card_W;
GGIW_.b = GGIW(end-f).b + 1;

z_bar = mean(W,2);
epsilon = z_bar - model.measmodel.h(GGIW(end-f).m);
H = model.measmodel.H(GGIW(end-f).m);

X_hat = GGIW(end-f).V/(GGIW(end-f).v - 2*d - 2);
X_hat = (X_hat + X_hat')/2;

S = H*GGIW(end-f).P*H' + X_hat/card_W;
S = (S + S')/2;

Vs= chol(S); det_S= prod(diag(Vs))^2; inv_sqrt_S= inv(Vs); iS= inv_sqrt_S*inv_sqrt_S';

K = GGIW(end-f).P*H'*iS;

GGIW_.m = GGIW(end-f).m + K*epsilon;
GGIW_.P = GGIW(end-f).P - K*H*GGIW(end-f).P;

temp = (W - z_bar);
Z = temp*temp';

X_sqrt = sqrtm_2by2(X_hat);
S_sqrt_inv = sqrtm(iS);
N = X_sqrt*S_sqrt_inv*(epsilon*epsilon')*S_sqrt_inv'*X_sqrt';

GGIW_.v = GGIW(end-f).v + card_W;
GGIW_.V = GGIW(end-f).V + N + Z;
GGIW_.model = length(model.motionmodel.F)-f;

%

lik = (GGIW(end-f).v-d-1)/2*log(det2(GGIW(end-f).V)) - (GGIW_.v-d-1)/2*log(det2(GGIW_.V))...
        + gamma2ln((GGIW_.v-d-1)/2) - gamma2ln((GGIW(end-f).v-d-1)/2)...
        + log(det2(X_hat))/2 - log(det_S)/2 + gammaln(GGIW_.a)...
        - gammaln(GGIW(end-f).a) + GGIW(end-f).a*log(GGIW(end-f).b) - GGIW_.a*log(GGIW_.b)...
        - (card_W*log(pi)+log(card_W))*d/2;
    
GGIW(end-f) = GGIW_;

end

