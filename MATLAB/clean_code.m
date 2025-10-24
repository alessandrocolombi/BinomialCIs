n_vec=[[500:100:2000] [2000:500:10000]]';
m=10000;
p=rand(m,1);
alph=0.05;
benchmark=log(m/alph)./n_vec;
mine_analytical=zeros(length(n_vec),1);
lb=zeros(length(n_vec),1);
ok_vec=zeros(length(n_vec),1);
num_of_exp=100;
for n_ind=1:length(n_vec)
    n=n_vec(n_ind)
    for exp_ind=1:num_of_exp
        r = binornd(n,p)/n;
        current_lb=0;
        if sum(r==0)>0
            current_lb=log(sum(r==0)/alph)/n;
        end
        lb(n_ind)=lb(n_ind)+current_lb;

        a1=0.99*alph;
        a2=0.01*alph;
        a=(log(n))^2;
        m_p=sum((1-r).^a)+a*sqrt(m*log(1/a2)/n);
        current_analytical=(log(m_p/a1)/(n-a));
        mine_analytical(n_ind)=mine_analytical(n_ind)+current_analytical;
        
        if sum(r==0)>0
            ok=0;
            if max(p(r==0))<=current_analytical
                ok=1;
            end
            ok_vec(n_ind)=ok_vec(n_ind)+ok;
        else
            ok_vec(n_ind)=ok_vec(n_ind)+1;
        end
    end

end
ok_vec=ok_vec/num_of_exp;
figure(2)
hold on
plot(n_vec,benchmark,'r')
plot(n_vec,mine_analytical/num_of_exp,'k')
plot(n_vec,lb/num_of_exp,'b--')

hold off
xlabel('n')
ylabel('bound')
%legend('benchmark','ours analytical')
legend('benchmark','ours analytical','lb')
box on
bla=1;