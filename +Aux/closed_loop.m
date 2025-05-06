function [x,u] = closed_loop(K,L,xo,Hor,Fig_name)

    x(:,1) = xo;
    
    for i = 1:Hor
        x(:,i+1) = L(:,:,Hor+1-i)*x(:,i);
        u(:,i) = K(:,:,Hor+1-i)*x(:,i);
    end

    figure('Name',Fig_name,'Position',[750 0 522 468])
    plot(1:Hor+1,x(1,:),'--k','LineWidth',2)
    hold on
    plot(1:Hor+1,x(2,:),':','Color',[0.5,0.5,0.5],'LineWidth',2)
    plot(1:Hor+1,x(3,:),'-.k','LineWidth',2)
    xlabel('i','fontsize',18,'fontweight','b')
    ylabel('x_{i}','fontsize',18,'fontweight','b')
    title('Closed-Loop System','fontsize',18)
    legend('x_{1,i}','x_{2,i}','x_{3,i}')
    set(gca,'FontSize',18,'FontWeight','bold')
    axis([1 Hor min(min(x))-0.1*abs(min(min(x))) ...
    max(max(x))+0.1*abs(max(max(x)))])
    grid on

    figure('Name',Fig_name,'Position',[750 600 522 468])
    plot(1:Hor,u(1,:),'-.k','LineWidth',2)
    hold on
    plot(1:Hor,u(2,:),':','Color',[0.5,0.5,0.5],'LineWidth',2)
    xlabel('i','fontsize',18,'fontweight','b')
    ylabel('u_{i}','fontsize',18,'fontweight','b')
    title('Control Action','fontsize',18)
    legend('u_{1,i}','u_{2,i}')
    set(gca,'FontSize',18,'FontWeight','bold')
    axis([1 Hor min(min(u))-0.1*abs(min(min(u))) ...
    max(max(u))+0.1*abs(max(max(u)))])
    grid on
end