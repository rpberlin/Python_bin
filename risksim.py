import random
import numpy as np
import matplotlib.pyplot as plt

def get_distrib_pt( x, cdf, cdf_val):
    if len(x) != len(cdf):
        return 0
    if cdf_val <= cdf[0]:
        return x[0]
    elif cdf_val >= cdf[-1]:
        return x[-1]
    for i in range(1,len(cdf)):
        x0 = cdf[i-1]
        x1 = cdf[i] 
        if cdf_val >= x0 and cdf_val < x1:
            y0 = x[i-1]
            y1 = x[i]
            return y0 + (y1-y0)*(cdf_val - x0)/(x1-x0)
    return 0
        
        

def dicerolls(n):
    res = [random.randint(1,6) for i in range(n)]
    res.sort()
    res.reverse()
    return res

def sim_blitz(a,d):
    #print(a,d)
    while a > 0 and d > 0:
        arisk = min(3,a)
        drisk = min(d,2)
        arolls = dicerolls(arisk)
        drolls = dicerolls(drisk)
        aloss = 0
        dloss = 0 
        n_battles = min(arisk,drisk)
        for i in range(0,n_battles):
            if arolls[i] <= drolls[i]:
                aloss += 1
            else:
                dloss += 1
        a -= aloss
        d -= dloss
        #print(a,d)
    return a, d
    

if __name__ == '__main__':
    n_reps = 1
    n_trials = 20000
    n_att_list = range(1,15)
    n_def = 5
    cdf_low = 0.2
    cdf_high = 0.8
    x_low_list = []
    x_mid_list = []
    x_high_list = []


    for n_att in n_att_list:
        print('N_att: ',n_att)
        xc_list = range(0,n_att+2)
        for j in range(n_reps):
            p_list = [0 for i in range(n_att+1)]
            for i in range(n_trials):
                a_rem, d_rem = sim_blitz(n_att,n_def)
                p_list[a_rem] += 1/n_trials
            c_list = [np.sum(p_list[0:i+1]) for i in range(0,len(p_list))]
            c_list.append(1.0)
            x_low = get_distrib_pt(xc_list, c_list, cdf_low)
            x_mid = get_distrib_pt(xc_list, c_list, 0.5)
            x_high = get_distrib_pt(xc_list,c_list,cdf_high)

            x_low_list.append(x_low)
            x_mid_list.append(x_mid)
            x_high_list.append(x_high)
            #plt.plot(range(0,n_att+2),c_list)
            #plt.plot([0,n_att+1],[cdf_low,cdf_low])
            #plt.plot([0,n_att+1],[cdf_high,cdf_high])
            #plt.plot([x_low,x_mid,x_high],[cdf_low,0.5,cdf_high],'o')
    plt.plot(n_att_list,x_low_list,label=f'CDF = {cdf_low}')
    plt.plot(n_att_list,x_mid_list,label='mid')
    plt.plot(n_att_list,x_high_list,label=f'CDF = {cdf_high}')
    plt.xlabel('Attacking Armies (-)')
    plt.ylabel('Remaining Armies of Attacker (-)')
    plt.title(f'Armies Remaining After Attacking {n_def} Defenders')
    plt.legend()
    plt.show()

    dr_da_low =[0]
    dr_da_mid =[0]
    dr_da_high =[0]
    for i in range(1,len(n_att_list)):
        dr_da_low.append(x_low_list[i]-x_low_list[i-1])
        dr_da_mid.append(x_mid_list[i]-x_mid_list[i-1])
        dr_da_high.append(x_high_list[i]-x_high_list[i-1])
    plt.plot(n_att_list,dr_da_low,label='low')
    plt.plot(n_att_list,dr_da_mid,label='mid')
    plt.plot(n_att_list,dr_da_high,label='high')
    plt.legend()
    plt.show()
    print('Done')

    
