function isOcclusion=isOcclusion(Line)
    isOcclusion=zeros(1,length(Line));
    for l=1:length(Line) %检验遮挡
        CenterPoint0=(Line{l}.start+Line{l}.end)/2;
        k0=(CenterPoint0(2)/CenterPoint0(1));%待检验线段中点与雷达连线的斜率
        for m=1:length(Line)%判断是否被其他线段遮挡
            if l~=m
                CenterPoint1=(Line{m}.start+Line{m}.end)/2;
                k1=Line{m}.k;
                x1=(CenterPoint1(2)-k1*CenterPoint1(1))/(k0-k1);
                if x1>=min(Line{m}.start(1),Line{m}.end(1)) &&  ...
                        x1<=max(Line{m}.start(1),Line{m}.end(1))...
                        && x1<CenterPoint0(1)
                    isOcclusion(l)=1;%标记线条
                    break;
                end
            end
        end
    end
end