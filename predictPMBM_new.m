function [PPP,MBM] = predictPMBM_new(PPP,MBM,model)
PPP.w = PPP.w + log(model.Ps);
PPP.GGIW = arrayfun(@(x) predictGGIWPPP(x,model), PPP.GGIW);

PPP.w = [PPP.w;log(model.birth.w)];
PPP.GGIW = [PPP.GGIW;model.birth.GGIW];


% Predict MBM
n_track = length(MBM.track);
for i = 1:n_track
    nh = length(MBM.track{i});
    for h = 1:nh
        if MBM.track{i}(h).Bern.w_death(end) >= model.threshold_s
            MBM.track{i}(h).Bern.GGIW = predictGGIW(MBM.track{i}(h).Bern.GGIW,model);
            MBM.track{i}(h).Bern.t_death = [MBM.track{i}(h).Bern.t_death MBM.track{i}(h).Bern.t_death(end)+1];
            MBM.track{i}(h).Bern.w_death = [MBM.track{i}(h).Bern.w_death(1:end-1) MBM.track{i}(h).Bern.w_death(end)...
                *(1-model.Ps) MBM.track{i}(h).Bern.w_death(end)*model.Ps];
            
        end
    end
end
%Forecast extension component
for i = 1:length(MBM.LineSegmentSet)
    nh = length(MBM.LineSegmentSet{i});
    for h = 1:nh
        for j=1:length(MBM.LineSegmentSet{i}(h).LineSegment)
            Bern=MBM.LineSegmentSet{i}(h).LineSegment{j}.Bern;
            Bern.GGIW = predictGGIW(Bern.GGIW,model);
            Bern.t_death = [Bern.t_death Bern.t_death(end)+1];
            Bern.w_death = [Bern.w_death(1:end-1) Bern.w_death(end)...
                *(1-model.Ps) Bern.w_death(end)*model.Ps];
            MBM.LineSegmentSet{i}(h).LineSegment{j}.Bern=Bern;
        end
    end
end
