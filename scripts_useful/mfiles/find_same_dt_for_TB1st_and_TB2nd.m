clear all
close all

t_dim=7;
your_desired_value=200;
search_diapazon=( round(your_desired_value-your_desired_value/2) : your_desired_value+1000);
% search_diapazon=2460:2550;
search_diapazon=500:600;
suitable_nts= {};

for i=search_diapazon
    suitable_nt_1st=[];
    suitable_nt_2nd=[];
    nt_2nd=i;
    nt_1st=2*i;
    remain_1st=remain(nt_1st,t_dim);
    remain_2nd=remain(nt_2nd,t_dim);
%     disp([remain_1st,remain_2nd])
    if (remain_1st==0) && (remain_2nd==0)
       suitable_nts=[suitable_nts,[nt_1st/2,nt_2nd]];
    else
        if (remain_1st==0)
            suitable_nt_1st=nt_1st/2;
            suitable_nt_2nd=nt_2nd+(t_dim+1)*2-remain_2nd;
        
        elseif (remain_2nd==0)
            suitable_nt_2nd=nt_2nd;
            suitable_nt_1st=(nt_1st+(t_dim+1)*2-remain_1st)/2;
        else
            suitable_nt_1st=(nt_1st+(t_dim+1)*2-remain_1st)/2;
            suitable_nt_2nd=nt_2nd+(t_dim+1)*2-remain_2nd;
        end
        disp([i,suitable_nt_1st,suitable_nt_2nd])
        if suitable_nt_1st==suitable_nt_2nd
            disp('!!!!!!')
            disp([i,suitable_nt_1st,suitable_nt_2nd])
            suitable_nts=[suitable_nts,[suitable_nt_1st,suitable_nt_2nd]];
        end
    end
end
ss=1
r_1st = remain(206,t_dim)
r_2nd = remain(206,t_dim)
% nt2 = p->nt + (p->t_dim+1)*2 - remain;
% int remain = (p->nt-2)%((p->t_dim+1)*2);

function r = remain(nt,t_dim)
    r=rem(nt-2,(t_dim+1)*2);
end
