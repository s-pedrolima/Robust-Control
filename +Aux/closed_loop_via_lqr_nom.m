function [x_mean,u_mean,J_mean] = closed_loop_via_lqr_nom(F, E_F, G, E_G, H, Q, R, K, xo, Hor)

    for j = 1:100
        x(:,1,j) = xo;
        J(:,1,j) = 0;
        
        for i = 1:Hor
            Delta = -1 + 2*rand(1);
            Delta_F = H*Delta*E_F;
            Delta_G = H*Delta*E_G;
            u(:,i,j) = -K*x(:,i,j);
            x(:,i+1,j) = ((F+Delta_F) - (G+Delta_G)*K)*x(:,i,j);
            J(:,i+1,j) = J(:,i,j) + ( x(:,i,j)'*Q*x(:,i,j) + ...
            u(:,i,j)'*R*u(:,i,j) );
        end
    end

    x_mean = mean(x,3);
    u_mean = mean(u,3);
    J_mean = mean(J,3);

    figure
    plot(1:Hor+1,x_mean(1,:),'--k','LineWidth',2)
    hold on
    plot(1:Hor+1,x_mean(2,:),':','Color',[0.5,0.5,0.5],'LineWidth',2)
    plot(1:Hor+1,x_mean(3,:),'-.k','LineWidth',2)
    xlabel('i','fontsize',18,'fontweight','b')
    ylabel('x_{i}','fontsize',18,'fontweight','b')
    title('Closed-Loop System (LQR)','fontsize',18)
    legend('x_{1,i}','x_{2,i}','x_{3,i}')
    set(gca,'FontSize',18,'FontWeight','bold')
    axis([1 Hor min(min(x_mean))-0.1*abs(min(min(x_mean))) ...
         max(max(x_mean))+0.1*abs(max(max(x_mean)))])
    grid on

    figure
    plot(1:Hor,u_mean(1,:),'-.k','LineWidth',2)
    hold on
    plot(1:Hor,u_mean(2,:),':','Color',[0.5,0.5,0.5],'LineWidth',2)
    xlabel('i','fontsize',18,'fontweight','b')
    ylabel('u_{i}','fontsize',18,'fontweight','b')
    title('Control Action (LQR)','fontsize',18)
    set(gca,'FontSize',18,'FontWeight','bold')
    legend('u_{1,i}','u_{2,i}')
    axis([1 Hor min(min(u_mean))-0.1*abs(min(min(u_mean))) ...
         max(max(u_mean))+0.1*abs(max(max(u_mean)))])
    grid on

    % figure(5)
    % plot(1:Hor+1,J_mean,'-.k','LineWidth',3)
    % xlabel('i','fontsize',14,'fontweight','b')
    % ylabel('J','fontsize',14,'fontweight','b')
    % title('Cost')
    % grid on
end