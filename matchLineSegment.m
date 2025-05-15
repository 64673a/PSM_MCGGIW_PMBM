function LineSegmentSet=matchLineSegment(LineSegmentSet,LineSegment,GGIW,model,time)
    if LineSegmentSet.detection==1
        l=model.lambda_l;
    else
        l=1;
    end
    [LineSegmentSet,~]=matchLineSegmentByLength(LineSegmentSet,LineSegment,GGIW,model,time,l);
end