function PseudoM=CreatePseudoMeasurement(LineSegment,model)
    precision=model.Ld_precision;
    M_num=0;
    for i=1:length(LineSegment)
        Astart=k2angle(LineSegment{i}.start(2),LineSegment{i}.start(1));
        Aend=k2angle(LineSegment{i}.end(2),LineSegment{i}.end(1));
        M_num=M_num+abs(Astart-Aend)/precision;
    end
    
    M_num=ceil(M_num/length(LineSegment));%-length(trackW(2,:));
    PseudoM=[];
    %M_num=;
    for i=1:length(LineSegment)
        LineW=[LineSegment{i}.start repmat(LineSegment{i}.end,[1,M_num])+[(LineSegment{i}.start-LineSegment{i}.end)/(M_num+1)]*[1:M_num] LineSegment{i}.end];
        PseudoM=[PseudoM LineW];
    end
end