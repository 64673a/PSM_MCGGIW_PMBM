function [unUseW,line]=getLineSegment(W,model,t)
line=[];
point=[];
unUseW=[];
nextPoint=1;%The number of points to try again after a failure
lambda_a=model.lambda_a;%Angle threshold

if ~isempty(W)
    %sort
    if length(W(1,:))>2
        measNum=size(W,2);
        for i=1:measNum
            W(3,i)=W(2,i)/W(1,i);%用角度排序
        end
        [~,W_index]=sort(W(3,:));%排序
        W=W(1:2,W_index);%保留量测
        
        %De-cluttering the measurement set
        W_l=zeros(1,measNum);
        for i=1:measNum
            if i==1
                W_l(i)=hypot(W(2,i+1)-W(2,i),W(1,i+1)-W(1,i));
            elseif i==measNum
                W_l(i)=hypot(W(2,i)-W(2,i-1),W(1,i)-W(1,i-1));
            else
                W_l(i)=min(hypot(W(2,i)-W(2,i-1),W(1,i)-W(1,i-1)),hypot(W(2,i+1)-W(2,i),W(1,i+1)-W(1,i)));
            end
        end
        if W_l(1)>2*mean(W_l(2:end))
            W=W(:,2:end);
            measNum=measNum-1;
        elseif W_l(end)>2*mean(W_l(1:end-1))
            W=W(:,1:end-1);
            measNum=measNum-1;
        end
            
        
        startPoint=1;
        lineIndex=1;
        endPoint=0;
        while endPoint~=measNum && startPoint<measNum%Cyclic group measurement
            endPoint=startPoint+1;
            point{lineIndex}=[W(:,startPoint) W(:,endPoint)];
            for toEndPoint=endPoint+1:measNum
                lineAk=k2angle(W(2,startPoint+1)-W(2,startPoint),W(1,startPoint+1)-W(1,startPoint));
                lineBk=k2angle((W(2,toEndPoint)-W(2,endPoint)),(W(1,toEndPoint)-W(1,endPoint)));
                if abs(lineAk-lineBk)<lambda_a || abs(lineAk-180)+lineBk<lambda_a || abs(lineBk-180)+lineAk<lambda_a
                    endPoint=toEndPoint;
                    point{lineIndex}=[point{lineIndex} W(:,endPoint)];
                else
                    anyNextPoint=0;
                    if toEndPoint < measNum-nextPoint+1 && toEndPoint~=measNum
                        for i=toEndPoint+1:toEndPoint+nextPoint
                            lineCk=k2angle((W(2,i)-W(2,endPoint)),(W(1,i)-W(1,endPoint)));
                            if abs(lineAk-lineCk)<lambda_a || abs(lineAk-180)+lineCk<lambda_a || abs(lineCk-180)+lineCk<lambda_a
                                endPoint=i;
                                point{lineIndex}=[point{lineIndex} W(:,endPoint)];
                                anyNextPoint=1;
                                break
                            end
                        end
                    end
                    if anyNextPoint==0
                        startPoint=endPoint;
                        lineIndex=lineIndex+1;
                        break;
                    end
                end
            end
        end
        
        %Delete the point set with only two points
        for i=1:length(point)
%             plot(point{i}(1,:),point{i}(2,:),'yo')
            if length(point{i})==2
                point{i}=[];
            end
        end
        point = point(cellfun(@(x) ~isempty(x), point));
        DeletelineTable=zeros(1,length(point));
        if length(point)>=2
            for i=1:length(point)

                if i==1
                    if point{i}(1,end)~=point{i+1}(1,1)
                        DeletelineTable(i)=1;
                    end
                end

                if i==length(point)
                    if point{i-1}(1,end)~=point{i}(1,1)
                        DeletelineTable(i)=1;
                    end
                end

                if i>1 && i<length(point)
                    if point{i-1}(1,end)~=point{i}(1,1) && point{i}(1,end)~=point{i+1}(1,1)
                        DeletelineTable(i)=1;
                    end
                end
            end
        end

        %TLS fit
        for i=1:length(point)
            n=length(point{i});
            X=point{i}(1,:)';
            Y=point{i}(2,:)';
            A=[ones(n,1) X];
            C=[A -Y];
            [~,R]=qr(C);
            C_R=R(2:end,2:end);
            [~,S,~]=svd(C_R);
            S=diag(S);
            p=-((A'*A-(S(2)^2)*[0 0;0 1])\(A')*(-Y))'*[0 1;1 0];
            dist=abs(p(1)*X-Y+p(2))./hypot(p(1),1);
            rms_error = sqrt(mean(dist.^2));
            dir = [1; p(1)]; dir = dir / norm(dir);
            center = mean([X Y])';
            V = [X Y]' - center;
            proj = dir' * V;
            proj_range = max(proj) - min(proj);

            lambda = 100;  % 惩罚项系数
            fit_quality = (n * proj_range) / (1 + lambda * rms_error^2);
            line{i}.p=p;
            line{i}.v=fit_quality;
            deviation(i)=fit_quality;
        end
        
        if ~isempty(line)
            if sum(DeletelineTable)==length(line) && ~isempty(line)%If all line segments are discrete, keep the two best line segments
                [~,idx]=sort(deviation);
                DeletelineTable(idx(end-1:end))=0;
            end
            for i=1:length(line)
                if DeletelineTable(i)==1
                    line{i}=[];
                    point{i}=[];
                end
            end
            line = line(cellfun(@(x) ~isempty(x), line));
            point = point(cellfun(@(x) ~isempty(x), point));
        end
        
        

            
        for i=1:length(line)
            if i==1
                m=point{i}(1,1);
                n=point{i}(2,1);
                a1=line{i}.p(1);
                b=-1;
                c1=line{i}.p(2);
                line{i}.start(1:2,1)=[(m+a1*n-a1*c1)/(a1^2+1),(a1^2*n+a1*m-b*c1)/(a1^2+1)];
                if length(line)>1
                    a2=line{i+1}.p(1);
                    c2=line{i+1}.p(2);
                    line{i}.end(1:2,1)=[(c2-c1)/(a1-a2),(c1*a2-c2*a1)/(a2-a1)];
                end
            end
            if i==length(line)
                m=point{i}(1,end);
                n=point{i}(2,end);
                a1=line{i}.p(1);
                c1=line{i}.p(2);
                line{i}.end(1:2,1)=[(m+a1*n-a1*c1)/(a1^2+1),(a1^2*n+a1*m+c1)/(a1^2+1)];
                if length(line)>1
                    a2=line{i-1}.p(1);
                    c2=line{i-1}.p(2);
                    line{i}.start(1:2,1)=[(c2-c1)/(a1-a2),(c1*a2-c2*a1)/(a2-a1)];
                end
            end
            if i~=1 && i~=length(line)
                line{i}.start=line{i-1}.end;
                a1=line{i}.p(1);
                c1=line{i}.p(2);
                a2=line{i+1}.p(1);
                c2=line{i+1}.p(2);
                line{i}.end(1:2,1)=[(c2-c1)/(a1-a2),(c1*a2-c2*a1)/(a2-a1)];
            end
        end
        
        %保证线段集连续Ensure that the segment set is continuous
        if length(line)==2 && point{1}(1,end)~=point{2}(1,1)
            line{3}=line{2};
            point{3}=point{2};
            
            m=point{1}(1,end);
            n=point{1}(2,end);
            a1=line{1}.p(1);
            c1=line{1}.p(2);
            line{1}.end(1:2,1)=[(m+a1*n-a1*c1)/(a1^2+1),(a1^2*n+a1*m+c1)/(a1^2+1)];
            line{2}.start=line{1}.end;
            
            m=point{3}(1,1);
            n=point{3}(2,1);
            a2=line{3}.p(1);
            c2=line{3}.p(2);
            line{3}.start(1:2,1)=[(m+a2*n-a2*c2)/(a2^2+1),(a2^2*n+a2*m+c2)/(a2^2+1)];
            line{2}.end=line{3}.start;
            
            
            line{2}.p(1)=(line{2}.end(2)-line{2}.start(2))/(line{2}.end(1)-line{2}.start(1));
            line{2}.v=(line{1}.v+line{3}.v)/2;
        end
        
        unUseW=W;
        for i=1:length(line)
            line{i}.matchCount=0;
            line{i}.birthTime=t;
            line{i}.k=line{i}.p(1);
            line{i}.point{t}=point{i};
            unUseW=setdiff(unUseW', point{i}','rows')';
            line{i}=rmfield(line{i},'p');
        end
    else
        unUseW=W;
    end
end
end