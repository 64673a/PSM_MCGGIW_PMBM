function LineW=getWforLine(Line,model)
    precision=model.Ld_precision;
    Astart=k2angle(Line.start(2),Line.start(1));
    Aend=k2angle(Line.end(2),Line.end(1));
    M_num=ceil(abs(Astart-Aend)/precision);
    LineW=[Line.start repmat(Line.end,[1,M_num])+[(Line.start-Line.end)/(M_num+1)]*[1:M_num] Line.end];
end