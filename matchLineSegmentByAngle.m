function [LineSegmentSet,sumCov]=matchLineSegmentByAngle(LineSegmentSet,LineSegment,GGIW,model,time,lambda_l)
    lambda_l=lambda_l*0.8;
    F=length(model.motionmodel.F);
    sumCov=inf;
    v=GGIW(end).m(3:4);
    if length(GGIW)>10
        v_=GGIW(end-5).m(3:4);
        c=5;
    else
        v_=GGIW(end-1).m(3:4);
        c=1;
    end
    theta=atan2d(v(1)*v_(2)-v(2)*v_(1),v(1)*v_(1)+v(2)*v_(2))/c;
    M=[cosd(theta),-sind(theta)
            sind(theta),cosd(theta)];
    v=(M*v);
    lambda_a=model.lambda_matchA;%角度差异阈值
    lambda_d=model.lambda_d;%长度差异阈值
    lmbda_m=1;%新生线段删除匹配次数阈值
    lambda_close=model.lambda_close;
if ~isempty(LineSegment)
    for i=1:length(LineSegmentSet.LineSegment) 
        matchLineSegment{i}.start=LineSegmentSet.LineSegment{i}.start+v;
        matchLineSegment{i}.end=LineSegmentSet.LineSegment{i}.end+v;
        matchLineSegment{i}.k=(matchLineSegment{i}.end(2)-matchLineSegment{i}.start(2))...
                /(matchLineSegment{i}.end(1)-matchLineSegment{i}.start(1));
    end
    isOcclusion=zeros(1,length(LineSegmentSet.LineSegment));
    matchLength=zeros(1,length(LineSegment));
    for l=1:length(matchLineSegment) %检验遮挡
        CenterPoint0=(matchLineSegment{l}.start+matchLineSegment{l}.end)/2;
        k0=(CenterPoint0(2)/CenterPoint0(1));%待检验线段中点与雷达连线的斜率
        for m=1:length(matchLineSegment)%判断是否被其他线段遮挡
            if l~=m
                CenterPoint1=(matchLineSegment{m}.start+matchLineSegment{m}.end)/2;
                k1=matchLineSegment{m}.k;
                x1=(CenterPoint1(2)-k1*CenterPoint1(1))/(k0-k1);
                if x1>=min(matchLineSegment{m}.start(1),matchLineSegment{m}.end(1)) &&  ...
                        x1<=max(matchLineSegment{m}.start(1),matchLineSegment{m}.end(1))...
                        && x1<CenterPoint0(1)
                    isOcclusion(l)=1;%标记线条
                    break;
                end
            end
        end
    end
    for i=1:length(LineSegment)
        lineAangle=k2angle(LineSegment{i});
        lineAlength=norm(LineSegment{i}.start-LineSegment{i}.end);
        for j=1:length(matchLineSegment)
            lineBangle=k2angle(matchLineSegment{j})-theta;
            lineBlength=norm(matchLineSegment{j}.start-matchLineSegment{j}.end);
            L=sum(abs((LineSegment{i}.start+LineSegment{i}.end)/2-(matchLineSegment{j}.start+matchLineSegment{j}.end)/2));
            lineD=min(lineAlength,lineBlength)/max(lineAlength,lineBlength);
            if (lineD>lambda_d && abs(lineBangle-lineAangle)<lambda_a && isOcclusion(j)~=1 && L<norm(v)*lambda_l) ||...
                (LineSegmentSet.Closed==1 && (j==1 || j==length(matchLineSegment)) && ...
                    lineD>lambda_d/2 && abs(lineBangle-lineAangle)<lambda_a*1.5 && isOcclusion(j)~=1 && L<norm(v)*lambda_l)
                if matchLength(1,i)==0
                    matchLength(1:2,i)=(LineSegment{i}.start+LineSegment{i}.end)/2-(matchLineSegment{j}.start+matchLineSegment{j}.end)/2;
                    matchLength(3,i)=sum(abs((LineSegment{i}.start+LineSegment{i}.end)/2-(matchLineSegment{j}.start+matchLineSegment{j}.end)/2));
                    matchLength(4,i)=i;
                    matchLength(5,i)=j;
                    matchLength(6,i)=abs(lineBangle-lineAangle);
                else
                    %匹配角度差异最小的两条线
                    if abs(lineBangle-lineAangle)<matchLength(6,i)
                        matchLength(1:2,i)=(LineSegment{i}.start+LineSegment{i}.end)/2-(matchLineSegment{j}.start+matchLineSegment{j}.end)/2;
                        matchLength(3,i)=sum(abs((LineSegment{i}.start+LineSegment{i}.end)/2-(matchLineSegment{j}.start+matchLineSegment{j}.end)/2));
                        matchLength(4,i)=i;
                        matchLength(5,i)=j;
                        matchLength(6,i)=abs(lineBangle-lineAangle);
                    end
                end
            end
        end
    end
    %删除为空的matchLength
    idx=0;
    for i=1:length(matchLength(1,:))
        idx=idx+1;
        if matchLength(1,idx)==0
            matchLength(:,idx)=[];
            idx=idx-1;
        end
    end
    if ~isempty(matchLength)    
        new_LineSegment=cell(0);
        matchLength=sortrows(matchLength',3)';
        dealt_L_sum=[0
                     0];
        dealt_a_sum=0;
        matchLineCount=0;
        matchEnd=0;
        v_sum=0;
        sumCov=inf;
        for i=find(matchLength(1,:))
            cov_xy=cov(matchLength(1,1:i),matchLength(2,1:i));
            cov_L=cov(matchLength(3,1:i));
            if abs(max(max(cov_xy)))>1 || cov_L>1
                break;
            else
                sumCov=mean(mean(cov_xy))+cov_L;
                matchLineCount=matchLineCount+1;
                matchEnd=i;
                v_sum=v_sum+LineSegment{matchLength(4,i)}.v;
            end
        end
        sumCov=sumCov/matchEnd+1;
        for i=1:matchEnd
            dealt_L=matchLength(1:2,i)*(LineSegment{matchLength(4,i)}.v/v_sum);
            dealt_a=(k2angle(LineSegment{matchLength(4,i)})-k2angle(matchLineSegment{matchLength(5,i)}))*(LineSegment{matchLength(4,i)}.v/v_sum);
            dealt_L_sum=dealt_L_sum+dealt_L;
            dealt_a_sum=dealt_a_sum+dealt_a;
        end
        dealt_a_sum=(dealt_a_sum-theta)/2;
        %删除多余matchLength值
        matchLength(:,matchEnd+1:end)=[];
        %匹配计数
        for i=1:length(LineSegmentSet.LineSegment)
            if ~isempty(find(matchLength(5,1:matchEnd)==i, 1))
                LineSegmentSet.LineSegment{i}.matchCount=LineSegmentSet.LineSegment{i}.matchCount+1;
            end
        end
        
        %移动set中的旧线段
        for i=1:length(LineSegmentSet.LineSegment)

                %位移
                new_LineSegment{end+1}.start=LineSegmentSet.LineSegment{i}.start+dealt_L_sum+v;
                new_LineSegment{end}.end=LineSegmentSet.LineSegment{i}.end+dealt_L_sum+v;
                %旋转
                x0=GGIW(end).m(1)+dealt_L_sum(1)+v(1);
                y0=GGIW(end).m(2)+dealt_L_sum(2)+v(2);
                    %start点
                w_x=new_LineSegment{end}.start(1);
                w_y=new_LineSegment{end}.start(2);
                w_x_turn=(w_x-x0)*cosd(dealt_a_sum)-(w_y-y0)*sind(dealt_a_sum)+x0;
                w_y_turn=(w_x-x0)*sind(dealt_a_sum)+(w_y-y0)*cosd(dealt_a_sum)+y0;
                new_LineSegment{end}.start(1)=w_x_turn;
                new_LineSegment{end}.start(2)=w_y_turn;
                    %end点
                w_x=new_LineSegment{end}.end(1);
                w_y=new_LineSegment{end}.end(2);
                w_x_turn=(w_x-x0)*cosd(dealt_a_sum)-(w_y-y0)*sind(dealt_a_sum)+x0;
                w_y_turn=(w_x-x0)*sind(dealt_a_sum)+(w_y-y0)*cosd(dealt_a_sum)+y0;
                new_LineSegment{end}.end(1)=w_x_turn;
                new_LineSegment{end}.end(2)=w_y_turn;
                %其它参数
                new_LineSegment{end}.k=(new_LineSegment{end}.start(2)-new_LineSegment{end}.end(2))/(new_LineSegment{end}.start(1)-new_LineSegment{end}.end(1));
                new_LineSegment{end}.matchCount=LineSegmentSet.LineSegment{i}.matchCount;
                new_LineSegment{end}.v=LineSegmentSet.LineSegment{i}.v;
                new_LineSegment{end}.birthTime=LineSegmentSet.LineSegment{i}.birthTime;
                new_LineSegment{end}.Bern=LineSegmentSet.LineSegment{i}.Bern;
                %记录链接点
                LinkPoints_new(:,i)=new_LineSegment{end}.start;
%                 figure(1);
%                 line([new_LineSegment{end}.start(1) new_LineSegment{end}.end(1)],[new_LineSegment{end}.start(2) new_LineSegment{end}.end(2)],'Color','green');
                %axis equal;
        end
        LinkPoints_new(:,end+1)=new_LineSegment{end}.end;
        
        %修正
        LinkPoints(:,1)=LineSegment{1}.start;
        for i=1:length(LineSegment)
            LinkPoints(:,i+1)=LineSegment{i}.end;
        end
        LinkPoints_new1=LinkPoints_new;
        LinkPoints_newTable=zeros(1,length(LinkPoints_new1(1,:)));
        PointMoveSum=zeros(2,length(LinkPoints_new1(1,:)));
        %判断是否要修正
        pointLock=zeros(1,length(LinkPoints_new1(1,:)));
        for i=1:length(isOcclusion)
            if isOcclusion(i)==1
                pointLock(i:i+1)=1;
            end
            if i>1 && i<length(isOcclusion)
                if isOcclusion(i+1)==0
                    pointLock(i)=0;
                elseif isOcclusion(i-1)==0
                    pointLock(i+1)=0;
                end
            end
        end
        for i=1:matchEnd
            v_sum=LineSegment{matchLength(4,i)}.v+LineSegmentSet.LineSegment{matchLength(5,i)}.v;
            if pointLock(matchLength(5,i))==0
                
                LinkPoints_newTable(matchLength(5,i))=1;
                
                
                PointMoveSum(:,matchLength(5,i))=PointMoveSum(:,matchLength(5,i))+LinkPoints_new1(:,matchLength(5,i))-...
                                              (LinkPoints_new(:,matchLength(5,i)) * (1-(LineSegment{matchLength(4,i)}.v/v_sum))...
                                             +LinkPoints(:,matchLength(4,i)) * (1-(LineSegmentSet.LineSegment{matchLength(5,i)}.v/v_sum)));
                                             %(LinkPoints_new(:,matchLength(5,i)) * 1/2 +LinkPoints(:,matchLength(4,i)) * 1/2);
                LinkPoints_new1(:,matchLength(5,i))=LinkPoints_new(:,matchLength(5,i)) * (1-(LineSegment{matchLength(4,i)}.v/v_sum))...
                                                    +LinkPoints(:,matchLength(4,i)) * (1-(LineSegmentSet.LineSegment{matchLength(5,i)}.v/v_sum));
                                             %(LinkPoints_new(:,matchLength(5,i)) * 1/2 +LinkPoints(:,matchLength(4,i)) * 1/2);
            end
            
            if pointLock(matchLength(5,i)+1)==0
                LinkPoints_newTable(matchLength(5,i)+1)=1;
                
                
                
                
                PointMoveSum(:,matchLength(5,i)+1)=PointMoveSum(:,matchLength(5,i)+1)+LinkPoints_new1(:,matchLength(5,i)+1)-...
                                                   (LinkPoints_new(:,matchLength(5,i)+1) * (1-(LineSegment{matchLength(4,i)}.v/v_sum))...
                                                   +LinkPoints(:,matchLength(4,i)+1) * (1-(LineSegmentSet.LineSegment{matchLength(5,i)}.v/v_sum)));
                
                                                    %(LinkPoints_new(:,matchLength(5,i)+1) * 1/2 +LinkPoints(:,matchLength(4,i)+1) * 1/2);
                
                
                
                LinkPoints_new1(:,matchLength(5,i)+1)=LinkPoints_new(:,matchLength(5,i)+1) * (1-(LineSegment{matchLength(4,i)}.v/v_sum))...
                                                      +LinkPoints(:,matchLength(4,i)+1) * (1-(LineSegmentSet.LineSegment{matchLength(5,i)}.v/v_sum));
                                                    %(LinkPoints_new(:,matchLength(5,i)+1) * 1/2 +LinkPoints(:,matchLength(4,i)+1) * 1/2);
            end
            %PointMoveSum
            new_LineSegment{matchLength(5,i)}.v=hypot(LineSegment{matchLength(4,i)}.v,LineSegmentSet.LineSegment{matchLength(5,i)}.v);
            %new_LineSegment{matchLength(5,i)}.v=(LineSegment{matchLength(4,i)}.v^2)/v_sum+LineSegmentSet.LineSegment{matchLength(5,i)}.v;
        end
        %未匹配线段修正
        if LinkPoints_newTable(1)==0
            for i=1:length(LinkPoints_newTable)
                if LinkPoints_newTable(i)==1
                    LinkPoints_new1(:,1:i-1)=LinkPoints_new1(:,1:i-1)-(PointMoveSum(:,i));
                    break
                end
            end
        end
        if LinkPoints_newTable(end)==0
            for i=length(LinkPoints_newTable):-1:1
                if LinkPoints_newTable(i)==1
                    LinkPoints_new1(:,i+1:end)=LinkPoints_new1(:,i+1:end)-(PointMoveSum(:,i));
                    break
                end
            end
        end
%         else
%            for i=1:length(LinkPoints_newTable)
%                 if LinkPoints_newTable==1
%                     LinkPoints_new1(:,i)
%                 end
%             end 
        %添加新的线段 

        headNewLine=0;
        tailNewLine=0;
        %newLineDeviation=[];
            %首部新线段
        if ~isempty(find(matchLength(5,1:matchEnd)==1, 1)) && LineSegmentSet.Closed==0%代表原set第一条线被匹配
            new_idx=matchLength(4,find(matchLength(5,1:matchEnd)==1, 1));
            if new_idx>1
                for i=1:new_idx-1
                    LinkPoints_new1=[LineSegment{new_idx-i}.start LinkPoints_new1];
                    headNewLine=headNewLine+1;
                end
            end
%         elseif LineSegmentSet.Closed==1
%             new_idx=matchLength(4,find(matchLength(5,1:matchEnd)==length(LineSegmentSet.LineSegment), 1));
%             if new_idx>1
%                 idx=1;
%                 for i=new_idx+1:length(LineSegment)
%                     LinkPoints_new1(:,end-idx)=LineSegment{new_idx-i}.start;
%                     idx=idx+1;
%                 end
%             end
        end

            %尾部新线段
        if ~isempty(find(matchLength(5,1:matchEnd)==length(LineSegmentSet.LineSegment), 1)) && LineSegmentSet.Closed==0%代表原set最后条线被匹配
            new_idx=matchLength(4,find(matchLength(5,1:matchEnd)==length(LineSegmentSet.LineSegment), 1));
            if new_idx<length(LineSegment)
                for i=new_idx+1:length(LineSegment)
                    LinkPoints_new1=[LinkPoints_new1 LineSegment{i}.end];
                    tailNewLine=tailNewLine+1;
                end
            end
%         elseif LineSegmentSet.Closed==1
%             new_idx=matchLength(4,find(matchLength(5,1:matchEnd)==length(LineSegmentSet.LineSegment), 1));
%             if new_idx<length(LineSegment)
%                 idx=0;
%                 for i=new_idx+1:length(LineSegment)
%                     LinkPoints_new1(:,2+idx)=LineSegment{i}.end;
%                     idx=idx+1;
%                 end
%             end
        end
        
        %赋值(头)
        new_LineSegment_head=[];
        for i=1:headNewLine
            new_LineSegment_head{i}.start=LinkPoints_new1(:,i);
            new_LineSegment_head{i}.end=LinkPoints_new1(:,i+1);
            new_LineSegment_head{i}.k=(new_LineSegment_head{i}.start(2)-new_LineSegment_head{i}.end(2))/...
                (new_LineSegment_head{i}.start(1)-new_LineSegment_head{i}.end(1));
            new_LineSegment_head{i}.matchCount=0;
            new_LineSegment_head{i}.birthTime=length(GGIW)+1;
            new_LineSegment_head{i}.v=LineSegment{i}.v;
            [Bern_L,lik_L]=getLineBern(new_LineSegment_head{i},model);
            Bern_L.t_birth = time;
            Bern_L.t_death = time;
            Bern_L.w_death = 1;
            new_LineSegment_head{i}.Bern=Bern_L;
            new_LineSegment_head{i}.lik=lik_L;
%             figure(1);
%             line([new_LineSegment_head{i}.start(1) new_LineSegment_head{i}.end(1)],[new_LineSegment_head{i}.start(2) new_LineSegment_head{i}.end(2)],'Color','blue','LineStyle','--');        
        end
        %赋值(旧)
        for i=1:length(new_LineSegment) 
            new_LineSegment{i}.start=LinkPoints_new1(:,headNewLine+i);
            new_LineSegment{i}.end=LinkPoints_new1(:,headNewLine+i+1);
            new_LineSegment{i}.k=(new_LineSegment{i}.start(2)-new_LineSegment{i}.end(2))/...
                (new_LineSegment{i}.start(1)-new_LineSegment{i}.end(1));
            Bern=new_LineSegment{i}.Bern;
            [Bern,lik]=detectionBern(Bern,getWforLine(new_LineSegment{i},model),model,length(model.motionmodel.F));
            new_LineSegment{i}.Bern=Bern;
            new_LineSegment{i}.lik=lik;
%             figure(1);
%             line([new_LineSegment{i}.start(1) new_LineSegment{i}.end(1)],[new_LineSegment{i}.start(2) new_LineSegment{i}.end(2)],'Color','black','LineStyle','--');
%             axis equal;
        end
        %赋值(尾)
        new_LineSegment_tail=[];
        idx=0;
        for i=tailNewLine:-1:1
            idx=idx+1;
            new_LineSegment_tail{idx}.start=LinkPoints_new1(:,end-i);
            new_LineSegment_tail{idx}.end=LinkPoints_new1(:,end-i+1);
            new_LineSegment_tail{idx}.k=(new_LineSegment_tail{idx}.start(2)-new_LineSegment_tail{idx}.end(2))/...
                (new_LineSegment_tail{idx}.start(1)-new_LineSegment_tail{idx}.end(1));
            new_LineSegment_tail{idx}.matchCount=0;
            new_LineSegment_tail{idx}.birthTime=length(GGIW)+1;
            new_LineSegment_tail{idx}.v=LineSegment{end-i+1}.v;
            [Bern_L,lik_L]=getLineBern(new_LineSegment_tail{idx},model);
            Bern_L.t_birth = time;
            Bern_L.t_death = time;
            Bern_L.w_death = 1;
            new_LineSegment_tail{idx}.Bern=Bern_L;
            new_LineSegment_tail{idx}.lik=lik_L;
            %             figure(1);
%             line([new_LineSegment_tail{idx}.start(1) new_LineSegment_tail{idx}.end(1)],[new_LineSegment_tail{idx}.start(2) new_LineSegment_tail{idx}.end(2)],'Color','yellow','LineStyle','--');
         end
        %合并
        LineSegmentSet.LineSegment=[new_LineSegment_head [new_LineSegment new_LineSegment_tail]];
        LineSegmentSet.Center=getLS_Center(LineSegmentSet.LineSegment);
        LineSegmentSet.allLineUnMatch=0;
        
        %删除线段
        deletLineTable=zeros(1,length(LineSegmentSet.LineSegment));
        for i=1:length(LineSegmentSet.LineSegment)
            if length(GGIW)+1-LineSegmentSet.LineSegment{i}.birthTime==2 && LineSegmentSet.LineSegment{i}.matchCount<=lmbda_m
                deletLineTable(i)=1;
                if i~=1 && i~=length(LineSegmentSet.LineSegment)
                    if i>length(LineSegmentSet.LineSegment)/2
                        deletLineTable(i+1:end)=1;
                    else
                        deletLineTable(1:i-1)=1;
                    end
                end
            end
        end
        for i=1:length(LineSegmentSet.LineSegment)
            if deletLineTable(i)==1 && sum(deletLineTable)~=length(LineSegmentSet.LineSegment)
                LineSegmentSet.LineSegment{i}=[];
            end
        end
        LineSegmentSet.LineSegment = LineSegmentSet.LineSegment(cellfun(@(x) ~isempty(x), LineSegmentSet.LineSegment));
        

        %%%%%%%收口
        if ~isempty(LineSegmentSet.LineSegment)
            perimeter=0;
            for i=1:length(LineSegmentSet.LineSegment)
                L=hypot(LineSegmentSet.LineSegment{i}.start(1)-LineSegmentSet.LineSegment{i}.end(1),...
                        LineSegmentSet.LineSegment{i}.start(2)-LineSegmentSet.LineSegment{i}.end(2));
                perimeter=perimeter+L;
            end
            closeLineL=hypot(LineSegmentSet.LineSegment{1}.start(1)-LineSegmentSet.LineSegment{end}.end(1),...
                        LineSegmentSet.LineSegment{1}.start(2)-LineSegmentSet.LineSegment{end}.end(2));
            perimeter=perimeter+closeLineL;%求出若封口后的周长,用于闭口判断
                    
            if closeLineL/perimeter<lambda_close
                closePoint=(LineSegmentSet.LineSegment{1}.start+LineSegmentSet.LineSegment{end}.end)/2;
                LineSegmentSet.LineSegment{1}.start=closePoint;
                LineSegmentSet.LineSegment{end}.end=closePoint;
                LineSegmentSet.LineSegment{1}.k=(LineSegmentSet.LineSegment{1}.start(2)-LineSegmentSet.LineSegment{1}.end(2))...
                    /(LineSegmentSet.LineSegment{1}.start(1)-LineSegmentSet.LineSegment{1}.end(1));
                LineSegmentSet.LineSegment{end}.k=(LineSegmentSet.LineSegment{end}.start(2)-LineSegmentSet.LineSegment{end}.end(2))...
                    /(LineSegmentSet.LineSegment{end}.start(1)-LineSegmentSet.LineSegment{end}.end(1));
                LineSegmentSet.LineSegment{end}.v=0;
                LineSegmentSet.LineSegment{1}.v=0;
                LineSegmentSet.Closed=1;
            else
                LineSegmentSet.Closed=0;
            end
        end
    else%如果所有线段都未匹配
        LineSegmentSet.allLineUnMatch=LineSegmentSet.allLineUnMatch+1;
        if LineSegmentSet.allLineUnMatch==3
            LineSegmentSet.LineSegment=LineSegment;
            for i=1:length(LineSegment)
                [Bern_L,lik_L]=getLineBern(LineSegment{i},model);
                Bern_L.t_birth = time;
                Bern_L.t_death = time;
                Bern_L.w_death = 1;
                LineSegmentSet.LineSegment{i}.Bern=Bern_L;
                LineSegmentSet.LineSegment{i}.lik=lik_L;
            end
            LineSegmentSet.Center=getLS_Center(LineSegment);
            LineSegmentSet.Closed=0;
            LineSegmentSet.allLineUnMatch=0;
        else
            for i=1:length(LineSegmentSet.LineSegment) 
                LineSegmentSet.LineSegment{i}.start=LineSegmentSet.LineSegment{i}.start+v;
                LineSegmentSet.LineSegment{i}.end=LineSegmentSet.LineSegment{i}.end+v;
                LineSegmentSet.LineSegment{i}.k=(LineSegmentSet.LineSegment{i}.end(2)-LineSegmentSet.LineSegment{i}.start(2))...
                        /(LineSegmentSet.LineSegment{i}.end(1)-LineSegmentSet.LineSegment{i}.start(1));
                [LineSegmentSet.LineSegment{i}.Bern,LineSegmentSet.LineSegment{i}.lik]=misdetectionBern(LineSegmentSet.LineSegment{i}.Bern,model,F);
            end
        end
    end
else
    for i=1:length(LineSegmentSet.LineSegment) 
        LineSegmentSet.LineSegment{i}.start=LineSegmentSet.LineSegment{i}.start+v;
        LineSegmentSet.LineSegment{i}.end=LineSegmentSet.LineSegment{i}.end+v;
        LineSegmentSet.LineSegment{i}.k=(LineSegmentSet.LineSegment{i}.end(2)-LineSegmentSet.LineSegment{i}.start(2))...
                /(LineSegmentSet.LineSegment{i}.end(1)-LineSegmentSet.LineSegment{i}.start(1));
        [LineSegmentSet.LineSegment{i}.Bern,LineSegmentSet.LineSegment{i}.lik]=misdetectionBern(LineSegmentSet.LineSegment{i}.Bern,model,F);
    end
end
deletLineTable=zeros(1,length(LineSegmentSet.LineSegment));
for i=1:length(LineSegmentSet.LineSegment)
    if length(GGIW)+1-LineSegmentSet.LineSegment{i}.birthTime==2 && LineSegmentSet.LineSegment{i}.matchCount<=lmbda_m
        deletLineTable(i)=1;
        if i~=1 && i~=length(LineSegmentSet.LineSegment) && i>length(LineSegmentSet.LineSegment)/2
            deletLineTable(i+1:end)=1;
        elseif i~=1 && i~=length(LineSegmentSet.LineSegment) && i<=length(LineSegmentSet.LineSegment)/2
            deletLineTable(1:i-1)=1;
        end
    end
end
for i=1:length(LineSegmentSet.LineSegment)
    if deletLineTable(i)==1 %&& sum(deletLineTable)~=length(LineSegmentSet.LineSegment)
        LineSegmentSet.LineSegment{i}=[];
    end
end
if ~isempty(LineSegmentSet.LineSegment)
    LineSegmentSet.LineSegment = LineSegmentSet.LineSegment(cellfun(@(x) ~isempty(x), LineSegmentSet.LineSegment));
end
end