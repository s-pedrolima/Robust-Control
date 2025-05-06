function [p,X,Y] = pos_polos(A,E_A,B,E_B,H,K,N,Fig_name,Titles)
    
    X = [];
    Y = [];
    
    for i = 1:N
        Delta = -1 + 2*rand(1);
        Delta_A = H*Delta*E_A;
        Delta_B = H*Delta*E_B;
        L = (A+Delta_A) + (B+Delta_B)*K;
        p(:,i) = eig(L);
        x = real(p(:,i));
        y = imag(p(:,i));
        X = [X ; x];
        Y = [Y ; y];
    end
    
    % Unit circle
    t = linspace(0,2*pi,N);
    h=0;
    k=0;
    r=1;
    coord_x = r*cos(t)+h;
    coord_y = r*sin(t)+k;
    
    % Plot
    figure('Name',Fig_name,'Position',[149 156 522 468])
    plot(coord_x,coord_y,'-.k','LineWidth',2);
    hold on
    plot(X,Y,'.k','MarkerSize',8)
    xlabel('Real','fontsize',12,'fontweight','b')
    ylabel('Imaginary','fontsize',12,'fontweight','b')
    title(Titles,'fontsize',14,'fontweight','b')
    axis equal;
    axis([-1.5 1.5 -1.5 1.5])
    set(gca,'FontSize',16,'FontWeight','bold')
    grid on

end