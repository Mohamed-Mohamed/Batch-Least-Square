function [ Gz ] = BatchLeastSquare ( u, y, na, nb, d, Ts, Plot )
% this function is used to estimate the parameter of the system without
% noise by useing BLS method
% estimated G(z)=z^d(b_0 z^nb+b_1 z^(nb-1)+...+b_nb)/(z^na+a_1 z^(na-1)+...+a_na)
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% inputs:
% u    : input required signal
% y      : output signal to the system
% na   : order of thr Den. of the transfer function of the system
% ab   : order of thr Num. of the transfer function of the system
% d     : order of delay of the system
% Ts     : sampling time
% Plot          : used to get the plot of system parameter estimation
%                   1 - if Plot = 0  --> no plot needed
%                   2 - if Plot = 1  -->  plot needed and Plot(2) == figuires number
%% outputs:
% Gz              : discreate transfer function

%% helper function
    function [ y_output ] = OutputEstimation ( A, B, de, U, Y, m )
        % order
        na1 = length(A) - 1;
        nb1 = length(B) -1;
        % coefficients
        B = B/A(1); % T ensure the first element (the highest power) is equal to 1
        A = A/A(1); % T ensure the first element (the highest power) is equal to 1
        A = -A(2:end); % to delete "1" at the first and -ve for the difference eqn
        % y(m)
        Y(m)=0;
        for j=1:na1
            if m-j >0
                Y(m)=Y(m)+A(j)*Y(m-j);
            end
        end
        for J=0:length(B)-1
            if m-J-de > 0
                Y(m)=Y(m)+B(J+1)*U(m-de-J);
            end
        end
        y_output=Y(m);
    end

%% Function body
for i=1:length(u)
    for j=1:na
        if i-j <=0
            phiT(i,j)=0;
        else
            phiT(i,j)=[-y(i-j)];
        end
    end
    for j=0:nb
        if i-j-d <= 0
            phiT(i,j+1+na)=0;
        else
            phiT(i,j+1+na)=[u(i-j-d)];
        end
    end
end
Theta_hat=inv(phiT'*phiT)*phiT'*y';
Gz=tf([Theta_hat(na+1:end)'],[1,Theta_hat(1:na)'],Ts);
%% plotting
if Plot(1)==1
    y1(1)=y(1);
    for k=1:length(y)-1
        [ y1_output ] = OutputEstimation ( [1,Theta_hat(1:na)'], [Theta_hat(na+1:end)'], d, u(1:k+1), y1(1:k), k+1 );
        y1(k+1)=y1_output;
    end
    figure(Plot(2));
    subplot(3,1,1:2)
    set(gcf,'color','w')
    plot((0:length(y)-1)*Ts,y,(0:length(y1)-1)*Ts,y1,'-o','linewidth',2);
    grid on;
    ylabel('y, y_e_s_t','fontsize',18);
    legend('y','y_e_s_t')
    subplot(3,1,3)
    plot((0:length(y1)-1)*Ts,abs(y-y1))
    grid on;
    xlabel('t(s)','fontsize',18);
    ylabel('y-y_e_s_t','fontsize',18);
    legend('error')
end
end