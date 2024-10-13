clear all
close all

dt=0.001;
t=(0:400)*dt;

% % ricker parameters
m_f0 = 8;
t0=1.5*sqrt(6.)/(pi*m_f0);
a  = pi*m_f0;
a2 = a  * a;
a4 = a2 * a2;
a6 = a4 * a2;

source_1st=zeros(length(t),1);
source_2nd=zeros(length(t),1);
for i=1:length(t)
    % 1st order
    deltaT=(t(i)-t0);
    deltaT2=deltaT  * deltaT;
    source_1st(i)=((1- 2* a2 * deltaT2)*exp(-a2*deltaT2))/dt;
    
    % 2nd order
    deltaT3=deltaT2 * deltaT;
    deltaT4=deltaT3 * deltaT;
    source_2nd(i)=(-6*a2*deltaT+4*a4*deltaT3)*exp(-a2*deltaT2);
end

font_sz=18;
figure;
plot(source_1st);
hold on
plot(source_2nd);
axis tight;
grid on
legend('source 1st','source 2nd');
set(gca,'FontSize',font_sz,'fontWeight','bold');
name='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/Apps/Overleaf/seg_23_alkhobar_workshop/Fig/source_1st_and_2nd_order.png'; delete(name)
print(name,'-dpng',['-r' num2str(400)]);

ss=1

figure
imagesc( squeeze(data(:,:,100)).' );
xlabel('X axis');
ylabel('Y axis');
colorbar


function r = remain(nt,t_dim)
    r=rem(nt-2,(t_dim+1)*2);
end
