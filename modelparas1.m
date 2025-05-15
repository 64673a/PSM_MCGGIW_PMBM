load('64673targetTracks.mat')    %从targetTracks.mat导入目标轨迹数据
load('64673.mat')         
%load('groundTruth_ManyTargets.mat')  
K = length(Scenario.Z{1});
model.K = K;

% Effective window length for the gamma prediction 伽马预测的有效窗口长度
w_e_gamma = 20;
% Effective window length for the extent prediction 扩展预测的有效窗口长度
w_e_extent = 10;

model.tao = 1/(log(w_e_extent)-log(w_e_extent-1));
model.eta = 1/(1-1/w_e_gamma);
model.Ts = 1;   %sampling interval 采样间隔
sigma_v = 0.1;  %standard deviation of motion noise 运动噪声的标准差
sigma_r = 0.1;  %standard deviation of measurement noise 测量噪声的标准差
model.motionmodel = motionmodel.immmodel(model.Ts,sigma_v);%调用motionmodel.m下的immmodel函数，传入采样间隔和标准差，返回d，F，Q，f 4个字段
                                                          %字段详细定义在motionmodel.m
model.measmodel = measmodel.cvmeasmodel(sigma_r);

%re-construct the data structure

% generate tracks (ground truth)
X = cell(K,1);
E = cell(K,1);
N = zeros(K,1);
groundTruth = cell(K,1); 
for targetnum = 1:length(targetTracks)          
    for k = targetTracks(targetnum).birthTime:targetTracks(targetnum).deathTime 
        targetstate = targetTracks(targetnum).x(1:model.motionmodel.d,k-targetTracks(targetnum).birthTime+1);
        targetextent = targetTracks(targetnum).X(:,:,k-targetTracks(targetnum).birthTime+1);
        X{k} = [X{k} targetstate];
        E{k} = cat(3,E{k},targetextent);
        N(k) = N(k) + 1;
    end
end

for k = 1:K
    groundTruth{k}.x = X{k};
    groundTruth{k}.X = E{k};
end

%target existence probability
model.Ps = 0.99;
%target detection probability
model.Pd = Scenario.detection_prob;

%range of the surveillance area
range_c = [-1 1;-1 1]*300;
%Poisson false alarm (clutter) rate 
lambda_c = Scenario.false_alarm_rate;
%Poisson clutter intensity
model.lambda_fa = lambda_c/prod(range_c(:,2)-range_c(:,1));

% target initial state
nbirths = 2;
xstart = zeros(model.motionmodel.d,nbirths);

xstart(:,1) = [500 250 0 0];
xstart(:,2) = [300 250 0 0];
%Birth model
d = 2;
model.birth.w = 0.5*ones(nbirths,1);
model.birth.GGIW = repmat(struct('a',250,'b',20,'m',[],'P',diag([5,5,5,5]),'v',10,'V',70*eye(2),'model',0),[nbirths,1]);
for i = 1:nbirths
    model.birth.GGIW(i).m = xstart(:,i);
end

% Gating parameters
Pg = 0.999;
model.gamma= chi2inv(Pg,model.measmodel.d);
model.Qd = 1 - model.Pd*Pg;
% Thresholds
model.threshold_r = 1e-2;   %existence probability of Bernoulli component
model.threshold_u = 1e-2;   %weight of mixture component in PPP 
model.threshold_w = 1e-2;   %1e-2 weight of global hypothesis (multi-Bernoulli) 
model.threshold_s = 1e-4;   %weight of the trajectory is still alive
model.recycle = 1e-1;       %recycling threshold 
model.merge = 4;            %merge threshold used to merge similar GGIWs 
model.M = 100;              %cap of number of MBM components in PMBM
model.num_iterations = 3;   %controls the number of iterations used in SO
model.max_repetition = 1;   %controls the number of iterations used in SO
model.Ld_precision = 0.06;  %Radar accuracy
model.lambda_a=45;          %Extraction angle threshold
model.lambda_matchA=30;     %Matching Angle Threshold
model.lambda_d=0.6;         %Match length threshold
model.lambda_close=0.1;     %Closure Threshold
model.lambda_l=0.7;         %Matching distance threshold
model.u_gauss=0.05;         %Gaussian noise intensity
%extract target state from Bernoulli components with existence probability
%larger than this threshold
model.exist_r = 0.5;        
