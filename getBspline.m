function [knot,ctrlpoints]=getBspline(LineSegment)
    ctrlpoints=[];
    knot=[0 0 1 1];
    for j=1:length(LineSegment)
        ctrlpoints=[ctrlpoints [LineSegment{j}.start LineSegment{j}.end]];
    end
    if length(ctrlpoints)>=3
        k_1=1;
        knot=[0 0 0 0];
        for k=1:length(ctrlpoints)-4
            knot=[knot k];
            k_1=k+1;
        end
        knot=[knot [k_1 k_1 k_1 k_1]];
    end
end