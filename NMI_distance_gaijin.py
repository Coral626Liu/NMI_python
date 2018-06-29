# coding=utf-8
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
import numpy as np
import math

def NMI(A,B):
    dict1=dict()
    dict2=dict()
    dict3=dict()
    first=DataFrame({"A":A,"B":B})
    tmp1=first["A"].groupby([first["A"],first["B"]])
    tmp2=first["A"].groupby([first["A"]])
    tmp3=first["B"].groupby([first["B"]])
    total=len(A)
    Hx=0
    Hy=0
    Ixy=0
    for i,j in tmp1:
        dict1[i]=float(len(j.values))/total
    for i,j in tmp2:        
        tmp=float(len(j.values))/total
        Hx += tmp*-1*math.log(tmp,2)
        dict2[i]=tmp
    for i,j in tmp3:
        tmp=float(len(j.values))/total
        Hy += tmp*-1*math.log(tmp,2)
        dict3[i]=tmp

    for i in dict1.keys():
        oppps=float(dict1[i])/(dict2[i[0]]*dict3[i[1]])
        Ixy+= dict1[i]*math.log(oppps,2)
    nmi=2*Ixy/(Hx+Hy)
    return nmi

def read_file_1(file,number1,number2):
    dat1=np.genfromtxt(file,skip_header=2, usecols=(1,6))
    df = DataFrame(dat1,columns=["frames","dist_%2d_%2d"%(number1,number2)])
    mins=df["dist_%2d_%2d"%(number1,number2)].groupby(df["frames"]).min() 
    return mins


def read_file_2(file):
    dat1=np.genfromtxt(file,skip_header=19, usecols=(0,1))                  
    df = DataFrame(dat1,columns=["frames","mindist"])
    return df
    #col=[]
#################################################################
#for i in xrange(10001):    
#    col.append("frame_%d"%i)

#mins_2=np.zeros(10001)
#for i in xrange(187,188):
#    for j in xrange(135,136):
#        mins=read_file("dist_%2d_%2d.data"%(i,j),i,j)
#        mins_2=np.vstack((mins_2,mins.values))   

#test=DataFrame(mins_2,columns=col)

#test.to_csv("./H_H6.csv")
#for i in xrange(10001):
 #  print test["frame_%d"%i].sort_values().index[0]
 #  print test["frame_%d"%i].sort_values().index[1]
 #   print test["frame_%d"%i].sort_values().index[2]
######################################################################
def get_auto_corr(timeSeries1_pre,timeSeries2_pre,k):
    """
     timeSeries is an array
    """
    l=len(timeSeries1_pre)
    timeSeries1=timeSeries1_pre[0:l-k]
    timeSeries2=timeSeries2_pre[k:]
    timeSeries1_mean=timeSeries1.mean()
    timeSeries2_mean=timeSeries2.mean()
    ###doubt
    timeSeries1_std= np.sqrt(timeSeries1_pre.var()*len(timeSeries1_pre))
    timeSeries2_std= np.sqrt(timeSeries2_pre.var()*len(timeSeries2_pre))
    auto_corr = 0
    for i in xrange(l-k):
        if timeSeries1_std == 0 or timeSeries2_std == 0:
            return 0
        else:
            tmp=(timeSeries1[i]-timeSeries1_mean)*(timeSeries2[i]-timeSeries2_mean)/(timeSeries1_std*timeSeries2_std)
            auto_corr = auto_corr + tmp
     
    return auto_corr
#####################################################################################
def plot_auto_corr(timeSeries1_pre,timeSeries2_pre,k,number1,number2):
    """
    k can not be beyound the length of timeSeries
    """
    timeSeriestimeSeries = pd.DataFrame(range(k))
    for i in xrange(1,k+1):
        timeSeriestimeSeries.loc[i-1] =get_auto_corr(timeSeries1_pre,timeSeries2_pre,i)
    plt.bar(range(1,len(timeSeriestimeSeries)+1),timeSeriestimeSeries[0].values)
    plt.savefig("./mind_hb_inter_%d_%d.png"%(number1,number2))
    plt.show()

mindist= read_file_2("../cap_ALA/mindist_10000_Ala.xvg")
mindist=mindist["mindist"]
for i in xrange(183,195):
    for j in xrange(130,182):
        mins=read_file_1("./dist_%2d_%2d.data"%(i,j),i,j)
        #tmp=df["hbnum_%d_%d"%(i,j)]
        nmi=NMI(mindist.values,mins.values)
        print i,j,"\t",nmi
############################################################################3
