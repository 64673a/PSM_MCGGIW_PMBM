function [GGIW] = predictGGIW(GGIW,model)

for i=1:length(model.motionmodel.F)
    d = 2;
    M = eye(2);
    e = exp(-model.Ts/model.tao);
    GGIW_{i}.v = 2*d + 2 + e*(GGIW(end).v - 2*d - 2);
    GGIW_{i}.V = e*M*GGIW(end).V*M';
    GGIW_{i}.a = GGIW(end).a/model.eta;
    GGIW_{i}.b = GGIW(end).b/model.eta;
    F = model.motionmodel.F{i}(GGIW(end).m);
    GGIW_{i}.m = F*GGIW(end).m;%状态预测步
    GGIW_{i}.P = F*GGIW(end).P*F' + model.motionmodel.Q;
    GGIW_{i}.model=i;
    if i~=1%旋转步
        if hypot(GGIW(end).m(3),GGIW(end).m(4))*hypot(GGIW_{i}.m(3),GGIW_{i}.m(4))~=0
            w=(GGIW(end).m(3)*GGIW_{i}.m(3)+GGIW(end).m(4)*GGIW_{i}.m(4))/(hypot(GGIW(end).m(3),GGIW(end).m(4))*hypot(GGIW_{i}.m(3),GGIW_{i}.m(4)));
            if GGIW_{i}.m(3)>GGIW(end).m(3)
                if GGIW_{i}.m(4)>0
                    rs=-acos(w);
                else
                    rs=acos(w);
                end
            else
                if GGIW_{i}.m(4)>0
                    rs=acos(w);
                else
                    rs=-acos(w);
                end
            end
            M=[cos(rs),-sin(rs)
                sin(rs),cos(rs)];
            GGIW_{i}.V = M*GGIW_{i}.V*M';
        end
    end
end

for i=1:length(model.motionmodel.F)
    GGIW = [GGIW;GGIW_{i}]; %存储预测值
end

end

