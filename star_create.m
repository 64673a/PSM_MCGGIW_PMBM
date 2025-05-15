
function star_create(numZ,detection_prob,false_alarm_rate,plotline)
rectangle=1;
star=1;
pcolor='c';
plinewidth=1.5;


precision=0.06;%Lidar accuracy
clutter=1;%clutter
maxTime=100;
Scenario.false_alarm_rate=false_alarm_rate; %Clutter rate
detection_prob_ld=0.99; %Single measurement detection probability
u_gauss=0.05; %Gaussian noise intensity
dplcount=0;

if star %星型目标
    targetTracks(1).x=[];
    targetTracks(1).X=[];
    for z=1:numZ
        detection_prob_flag=0;
        %Generate star convex lidar measurements 生成星凸型激光雷达量测
        center=[300 250]; %Center coordinates 中心坐标
        length=4; %Side length边长
        v=[8 0];%Initial x,y direction velocity 初始x,y方向速度
        

        p{1}=[center(1)-1.5*length,center(2)+0.5*length];
        p{2}=[center(1)-0.5*length,center(2)+0.5*length];
        p{3}=[center(1)-0.5*length,center(2)+1.5*length];
        p{4}=[center(1)+0.5*length,center(2)+1.5*length];
        p{5}=[center(1)+0.5*length,center(2)+0.5*length];
        p{6}=[center(1)+1.5*length,center(2)+0.5*length];
        p{7}=[center(1)+1.5*length,center(2)-0.5*length];
        p{8}=[center(1)+0.5*length,center(2)-0.5*length];
        p{9}=[center(1)+0.5*length,center(2)-1.5*length];
        p{10}=[center(1)-0.5*length,center(2)-1.5*length];
        p{11}=[center(1)-0.5*length,center(2)-0.5*length];
        p{12}=[center(1)-1.5*length,center(2)-0.5*length];
        p{13}=p{1};
        
        Scenario.Z{z}=cell(1,maxTime);
        
        
        for i=1:maxTime
            [1 z i]
            m=[center v]';
            
            if plotline
                figure(11);
                axis([100,700,100,700]);
                axis equal
                for j=2:13
                    plot([p{j-1}(1) p{j}(1)],[p{j-1}(2) p{j}(2)],'-','linewidth',plinewidth,'color',pcolor,'DisplayName','RealShape');
                    hold on;
                end
            end
            %存储真实数据
            if z == 1
                targetTracks(1).x(1:4,i)=m;
                pt=[];
                for j=1:12
                    t(1,:)=linspace(p{j}(1),p{j+1}(1),10000);
                    t(2,:)=linspace(p{j}(2),p{j+1}(2),10000);
                    pt=[pt t];
                end
                X=cov(pt(1,:),pt(2,:));
                targetTracks(1).X(:,:,i)=X;
            end
            
            if rand()<detection_prob || detection_prob_flag==1 || i<=3
                detection_prob_flag=0;
                k_=atand(center(2)/center(1));
                for j=k_-3:precision:k_+3 %激光雷达扫描扇形的范围
                    k1=tand(j);
                    temporary=[];
                    for k=2:13
                        if((p{k}(1)-p{k-1}(1))~=0)
                            k2=(p{k}(2)-p{k-1}(2))/(p{k}(1)-p{k-1}(1));
                            u=[(-p{k}(1)*k2+p{k}(2))/(k1-k2),k1*(((-p{k}(1)*k2+p{k}(2))/(k1-k2)))]';
                            u(3)=sqrt(u(1)^2+u(2)^2);
                            if (u(1)<=max(p{k}(1),p{k-1}(1))+0.01 && u(1)>=min(p{k}(1),p{k-1}(1))-0.01 ...
                                    && u(2)<=max(p{k}(2),p{k-1}(2))+0.01 && u(2)>=min(p{k}(2),p{k-1}(2))-0.01)
                                temporary=[temporary,u];
                            end
                        end
                        if((p{k}(1)-p{k-1}(1))==0)
                            u=[p{k}(1),k1*p{k}(1)]';
                            u(3)=sqrt(u(1)^2+u(2)^2);
                            if (u(1)<=max(p{k}(1),p{k-1}(1))+0.01 && u(1)>=min(p{k}(1),p{k-1}(1))-0.01 ...
                                    && u(2)<=max(p{k}(2),p{k-1}(2))+0.01 && u(2)>=min(p{k}(2),p{k-1}(2))-0.01)
                                temporary=[temporary,u];
                            end
                        end
                    end
                    dpl=rand();
                    if dpl>=detection_prob_ld && ~isempty(temporary)
                        dplcount=dplcount+1;
                    end
                    if ~isempty(temporary) && dpl<=detection_prob_ld
                        temporarys=temporary';
                        temporarys=sortrows(temporarys,3);
                        %hold on;
                        temporarys(1,1)=temporarys(1,1)+randn()*u_gauss;
                        temporarys(1,2)=temporarys(1,2)+randn()*u_gauss;
                        Scenario.Z{z}{i}=[Scenario.Z{z}{i},temporarys(1,1:2)'];
                    end
                end
            else
                detection_prob_flag=1;
            end
            %转态改变
            center=center+v;
            rs=4;
            if i<=20
                rs=3;
            end

            for k=1:13
                p{k}=p{k}+v;
                rx=(p{k}(1)-center(1))*cosd(rs)-(p{k}(2)-center(2))*sind(rs)+center(1);
                ry=(p{k}(1)-center(1))*sind(rs)+(p{k}(2)-center(2))*cosd(rs)+center(2);
                p{k}(1)=rx;
                p{k}(2)=ry;
            end
            %转动
            
            M=[cosd(rs),-sind(rs)
                sind(rs),cosd(rs)];%转动矩阵
            v=(M*v')';%速度转动
            if clutter
                clutterNum=poissrnd(Scenario.false_alarm_rate);%杂波生成
                for j=1:clutterNum
                    cx=100+600*rand();
                    cy=100+600*rand();
                    Scenario.Z{z}{i}=[Scenario.Z{z}{i},[cx,cy]'];
                end
            end
        end
        
        targetTracks(1).birthTime=1;
        targetTracks(1).deathTime=100;
    end
end


if rectangle
    for z=1:numZ
        detection_prob_flag=0;
        T=1;%扫描间隔
        maxTime=100;%总时间
        s=-3;%每次逆时针旋转随机0~s度
        F = [1, T, 0, 0
            0, 1, 0, 0
            0, 0, 1, T
            0, 0, 0, 1]; %目标转移矩阵
        W=8;
        H=6;
        V=[-8 0];%初始x方向速度
        a=[500-(W/2),V(1),250-(H/2),V(2)]';%分别为[x轴起始坐标，x方向速度，y起始坐标，y方向速度] 矩形左下角坐标及参数
        b=[a(1)+W,a(2),a(3),a(4)]';
        c=[a(1)+W,a(2),a(3)+H,a(4)]';
        d=[a(1),a(2),a(3)+H,a(4)]';
        
        for i=1:T:maxTime
            [2 z i]
            %绘图
            if plotline
                figure(11);
                axis([100,700,100,700]);

                plot([a(1) b(1)],[a(3) b(3)],'-','linewidth',plinewidth,'color',pcolor,'DisplayName','RealShape');
                plot([b(1) c(1)],[b(3) c(3)],'-','linewidth',plinewidth,'color',pcolor,'DisplayName','RealShape');
                plot([c(1) d(1)],[c(3) d(3)],'-','linewidth',plinewidth,'color',pcolor,'DisplayName','RealShape');
                plot([d(1) a(1)],[d(3) a(3)],'-','linewidth',plinewidth,'color',pcolor,'DisplayName','RealShape');
                hold on;
            end
            p1(1,:)=[a(1),b(1),c(1),d(1),a(1)];
            p1(2,:)=[a(3),b(3),c(3),d(3),a(3)];
            
            if z==1
                targetTracks(2).x(1:4,i)=[(a(1)+c(1))/2,(a(3)+c(3))/2,a(2),a(4)]';
                pt=[];
                for j=1:4
                    t(1,:)=linspace(p1(1,j),p1(1,j+1),10000);
                    t(2,:)=linspace(p1(2,j),p1(2,j+1),10000);
                    pt=[pt t];
                end
                X=cov(pt(1,:),pt(2,:));
                targetTracks(2).X(:,:,i)=X;
            end
            
            
            if rand()<detection_prob || detection_prob_flag==1 || i<=3
                detection_prob_flag=0;
                k_=atand((a(3)+c(3))/(a(1)+c(1)));
                for j=k_-5:precision:k_+5 %激光雷达扫描扇形的范围
                    k1=tand(j);
                    xs=[a(1),b(1),c(1),d(1),a(1)];
                    ys=[a(3),b(3),c(3),d(3),a(3)];
                    temporary=[];
                    for k=2:5
                        if((xs(k)-xs(k-1))~=0)
                            k2=(ys(k)-ys(k-1))/(xs(k)-xs(k-1));
                            u=[(-xs(k)*k2+ys(k))/(k1-k2),k1*(((-xs(k)*k2+ys(k))/(k1-k2)))]';
                            u(3)=sqrt(u(1)^2+u(2)^2);
                            if (u(1)<=max(xs(k),xs(k-1))+0.01 && u(1)>=min(xs(k),xs(k-1))-0.01 && u(2)<=max(ys(k),ys(k-1))+0.01 && u(2)>=min(ys(k),ys(k-1))-0.01)
                                temporary=[temporary,u];
                            end
                        end
                        if((xs(k)-xs(k-1))==0)
                            u=[xs(k),k1*xs(k)]';
                            u(3)=sqrt(u(1)^2+u(2)^2);
                            if (u(1)<=max(xs(k),xs(k-1))+0.01 && u(1)>=min(xs(k),xs(k-1))-0.01 && u(2)<=max(ys(k),ys(k-1))+0.01 && u(2)>=min(ys(k),ys(k-1))-0.01)
                                temporary=[temporary,u];
                            end
                        end
                    end
                    dpl=rand();
                    if dpl>=detection_prob_ld
                        dplcount=dplcount+1;
                    end
                    if ~isempty(temporary) && dpl<=detection_prob_ld
                        temporarys=temporary';
                        temporarys=sortrows(temporarys,3);
                        temporarys(1,1)=temporarys(1,1)+randn()*u_gauss;
                        temporarys(1,2)=temporarys(1,2)+randn()*u_gauss;
                        %hold on;
                        %                     plot(temporarys(1,1),temporarys(1,2),'*');
                        Scenario.Z{z}{i}=[Scenario.Z{z}{i},temporarys(1,1:2)'];
                    end
                end
            else
                detection_prob_flag=1;
            end
            %状态改变（主运动）
            a=F*a;
            b=F*b;
            c=F*c;
            d=F*d;
            %转动
            rx0=(a(1)+c(1))/2;
            ry0=(a(3)+c(3))/2;
            rs=-4;%+(rand()-0.5)*s*2;
            if i<=20
                rs=-3;
            end

            M=[cosd(rs),-sind(rs)
                sind(rs),cosd(rs)];%转动矩阵
            va=[a(2),a(4)]';
            va=M*va;
            a(2)=va(1);
            a(4)=va(2);
            b(2)=va(1);
            b(4)=va(2);
            c(2)=va(1);
            c(4)=va(2);
            d(2)=va(1);
            d(4)=va(2);
            ax=(a(1)-rx0)*cosd(rs)-(a(3)-ry0)*sind(rs)+rx0;
            a(3)=(a(1)-rx0)*sind(rs)+(a(3)-ry0)*cosd(rs)+ry0;
            bx=(b(1)-rx0)*cosd(rs)-(b(3)-ry0)*sind(rs)+rx0;
            b(3)=(b(1)-rx0)*sind(rs)+(b(3)-ry0)*cosd(rs)+ry0;
            cx=(c(1)-rx0)*cosd(rs)-(c(3)-ry0)*sind(rs)+rx0;
            c(3)=(c(1)-rx0)*sind(rs)+(c(3)-ry0)*cosd(rs)+ry0;
            dx=(d(1)-rx0)*cosd(rs)-(d(3)-ry0)*sind(rs)+rx0;
            d(3)=(d(1)-rx0)*sind(rs)+(d(3)-ry0)*cosd(rs)+ry0;
            a(1)=ax;
            b(1)=bx;
            c(1)=cx;
            d(1)=dx;
        end
        targetTracks(2).birthTime=1;
        targetTracks(2).deathTime=100;
    end
end













% figure(2)
% for i = 1:size(Scenario.Z{z},2)
%     %pause(0.1)
%     if ~isempty(Scenario.Z{z}{i})
%         plot(Scenario.Z{z}{i}(1,:),Scenario.Z{1}{i}(2,:),'b.');
%         hold on;
%     end
% end
disp('done')
Scenario.detection_prob=detection_prob;
save ('64673','Scenario') ;
save ('64673targetTracks','targetTracks');
end