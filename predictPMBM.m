function [PPP,MBM] = predictPMBM(PPP,MBM,model)

% Predict existing PPP 预测现有PPP
PPP.w = PPP.w + log(model.Ps);
PPP.GGIW = arrayfun(@(x) predictGGIWPPP(x,model), PPP.GGIW);

% Incorporate PPP birth 纳入PPP.birth
PPP.w = [PPP.w;log(model.birth.w)];
PPP.GGIW = [PPP.GGIW;model.birth.GGIW]; %因为默认有4个机场,所以在每次更新的时候把机场加回未检测目标内,方便新目标更新


% Predict MBM
n_track = length(MBM.track); %有n_track个track
for i = 1:n_track
    nh = length(MBM.track{i}); %每个track有n个局部假设
    for h = 1:nh
        if MBM.track{i}(h).Bern.w_death(end) >= model.threshold_s %轨迹存活权重大于预设的值
            MBM.track{i}(h).Bern.GGIW = predictGGIW(MBM.track{i}(h).Bern.GGIW,model); %预测该轨迹的GGIW
            MBM.track{i}(h).Bern.t_death = [MBM.track{i}(h).Bern.t_death MBM.track{i}(h).Bern.t_death(end)+1]; %当前时刻存在,死亡时刻+1s
            MBM.track{i}(h).Bern.w_death = [MBM.track{i}(h).Bern.w_death(1:end-1) MBM.track{i}(h).Bern.w_death(end)...
                *(1-model.Ps) MBM.track{i}(h).Bern.w_death(end)*model.Ps];%每个时刻目标的存活权重
            
        end
    end
end

end

