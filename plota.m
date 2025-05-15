clear;clc
close all;
load ('64673','measurement') ;
figure(2)
for i = 1:size(measurement,2)
    pause(0.1)
    if ~isempty(measurement{i})
        plot(measurement{i}(1,:),measurement{i}(2,:),'b.');
        axis([-100 300 -100 300]);
        hold on;
    end
end