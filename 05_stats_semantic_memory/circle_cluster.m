function circle_cluster(x_axis,y_axis,color,lineWidth);

%%%% this is the function to circle the significant clusters
%%%% x_axis are the corrdinate of all the x axis
%%%% y_axis are the corrdinate of all the y axis
%%%% color and lineWidth are the property of the line
%%%% by hui April-05-2012
for j = 1:length(x_axis)
    tmp_X = x_axis(j);
    tmp_Y = y_axis(j);
    tmp_1 = find(x_axis==(tmp_X + 1));
    if length(tmp_1)>0 && any(y_axis(tmp_1) == tmp_Y)
%         disp('Hi neibour!');
    else             
        x_cor = ones(1,10)*(tmp_X+0.5);%[(tmp_X+0.5),(tmp_X+0.5)];
        y_cor = linspace(tmp_Y-0.5,tmp_Y+0.5,10);%[(tmp_Y-0.5),(tmp_Y+0.5)];
        plot(x_cor,y_cor,color,'LineWidth',lineWidth);
        hold on
        clear x_cor y_cor
    end
    clear tmp_1
    tmp_1 = find(x_axis==(tmp_X - 1));
    if length(tmp_1)>0 && any(y_axis(tmp_1) == tmp_Y)
%         disp('Hi neibour!');
    else
        x_cor = ones(1,10)*(tmp_X-0.5);%[(tmp_X-0.5),(tmp_X-0.5)];
        y_cor = linspace(tmp_Y-0.5,tmp_Y+0.5,10);%[(tmp_Y-0.5),(tmp_Y+0.5)];
        plot(x_cor,y_cor,color,'LineWidth',lineWidth);
        hold on
        clear x_cor y_cor
    end
    clear tmp_1
    tmp_1 = find(x_axis==tmp_X);
    if length(tmp_1)>0 && any(y_axis(tmp_1) == (tmp_Y+1))
%         disp('Hi neibour!');
    else
        x_cor = linspace(tmp_X-0.5,tmp_X+0.5,10);%[(tmp_X-0.5),(tmp_X+0.5)];
        y_cor = ones(1,10)*(tmp_Y+0.5);%[(tmp_Y+0.5),(tmp_Y+0.5)];
        plot(x_cor,y_cor,color,'LineWidth',lineWidth);
        hold on
        clear x_cor y_cor
    end
    clear tmp_1
    tmp_1 = find(x_axis==tmp_X);
    if length(tmp_1)>0 && any(y_axis(tmp_1) == (tmp_Y-1))
%         disp('Hi neibour!');
    else
        x_cor = linspace(tmp_X-0.5,tmp_X+0.5,10);%[(tmp_X-0.5),(tmp_X+0.5)];
        y_cor = ones(1,10)*(tmp_Y-0.5);%[(tmp_Y-0.5),(tmp_Y-0.5)];
        plot(x_cor,y_cor,color,'LineWidth',lineWidth);
        hold on
        clear x_cor y_cor
    end
    clear tmp_X tmp_Y tmp_1
end