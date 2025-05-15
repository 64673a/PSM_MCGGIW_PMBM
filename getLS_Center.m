function LS_Center=getLS_Center(LineSegment)
    flag=0;
    Left=inf;
    Right=-inf;
    Up=-inf;
    Down=inf;
    for l=1:length(LineSegment) %检验遮挡
        if LineSegment{l}.start(1)<Left || LineSegment{l}.end(1)<Left
            Left=min(LineSegment{l}.start(1),LineSegment{l}.end(1));
        end
        if LineSegment{l}.start(1)>Right || LineSegment{l}.end(1)>Right
            Right=max(LineSegment{l}.start(1),LineSegment{l}.end(1));
        end
        if LineSegment{l}.start(2)>Up || LineSegment{l}.end(2)>Up
            Up=max(LineSegment{l}.start(2),LineSegment{l}.end(2));
        end
        if LineSegment{l}.start(2)<Down || LineSegment{l}.end(2)<Down
            Down=min(LineSegment{l}.start(2),LineSegment{l}.end(2));
        end
        
        CenterPoint0=(LineSegment{l}.start+LineSegment{l}.end)/2;
        k0=(CenterPoint0(2)/CenterPoint0(1));%待检验线段中点与雷达连线的斜率
        %line([0 CenterPoint0(1)],[0 CenterPoint0(2)]);
        for m=1:length(LineSegment)%判断是否被其他线段遮挡
            if l~=m &&  ~isempty(LineSegment{m})
                CenterPoint1=(LineSegment{m}.start+LineSegment{m}.end)/2;
                k1=LineSegment{m}.k;
                x1=(CenterPoint1(2)-k1*CenterPoint1(1))/(k0-k1);
                if x1>=min(LineSegment{m}.start(1),LineSegment{m}.end(1)) &&  ...
                        x1<=max(LineSegment{m}.start(1),LineSegment{m}.end(1))...
                        && x1<CenterPoint0(1)
                        flag=1;
                    break;
                end
            end
        end
    end
    if flag==1
        LS_Center(1)=(Left+Right)/2;
        LS_Center(2)=(Up+Down)/2;
    else
%         a=-(LineSegment{1}.start(2)-LineSegment{end}.end(2))/(LineSegment{1}.start(1)-LineSegment{end}.end(1));
%         b=1;
%         c=(-a)*LineSegment{1}.start(1)-LineSegment{1}.start(2);
%         for l=1:1:length(LineSegment)
%             x0=LineSegment{l}.start(1);
%             y0=LineSegment{l}.start(2);
%             LineSegment1{l}.start(1)=x0-2*a*((a*x0+b*y0+c)/(a^2+b^2));
%             LineSegment1{l}.start(2)=y0-2*b*((a*x0+b*y0+c)/(a^2+b^2));
%             x0=LineSegment{l}.end(1);
%             y0=LineSegment{l}.end(2);
%             LineSegment1{l}.end(1)=x0-2*a*((a*x0+b*y0+c)/(a^2+b^2));
%             LineSegment1{l}.end(2)=y0-2*b*((a*x0+b*y0+c)/(a^2+b^2));
%             
%             if LineSegment1{l}.start(1)<Left || LineSegment1{l}.end(1)<Left
%                 Left=min(LineSegment1{l}.start(1),LineSegment1{l}.end(1));
%             elseif LineSegment1{l}.start(1)>Right || LineSegment1{l}.end(1)>Right
%                 Right=max(LineSegment1{l}.start(1),LineSegment1{l}.end(1));
%             elseif LineSegment1{l}.start(2)>Up || LineSegment1{l}.end(2)>Up
%                 Up=min(LineSegment1{l}.start(2),LineSegment1{l}.end(2));
%             elseif LineSegment1{l}.start(2)<Down || LineSegment1{l}.end(2)<Down
%                 Down=max(LineSegment1{l}.start(2),LineSegment1{l}.end(2));
%             end
%         end
%         for l=1:length(LineSegment1)
%         line([LineSegment1{l}.start(1) LineSegment1{l}.end(1)],[LineSegment1{l}.start(2) LineSegment1{l}.end(2)]);
%         end
        LS_Center(1)=(Left+Right)/2;
        LS_Center(2)=(Up+Down)/2;
    end
%     plot(LS_Center(1),LS_Center(2),'r.');
end