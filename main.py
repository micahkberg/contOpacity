"""
Created on Tue Jan 19 13:42:41 2016

@author: iceberg

this is going to model continuous opacities of stars, someday
"""

#import stuff
import numpy as np
import matplotlib.pyplot as plt
plt.ion() #did i use this?
data = np.loadtxt("temp_pelog.dat")
correct=np.loadtxt("correct.dat")
import msvcrt as m  #see next comment
def wait():  #dont think i ever used this, was trying to get a plot to work inside of a function, but i think plotting from inside a function hangs the program somehow
    m.getch()
#define constants
sigma=5.67051e-5            #ergs s-1 cm-2 K-4, stefan-boltz
planck=6.6260755E-27        #erg seconds
boltz=1.380658E-16          #erg K-1
speed_light=2.99792458E10   #cm/sec
pmass=1.67262158E-24        #grams
e    =2.718281828459        #eulers i
chi     =2.179E-11          #hydrogen ionization in ergs
m_H     =1.6735575E-24      #mass of hydrogen in grams
sigma_t =6.655E-25          #cm^2
B       =0.1                #n_He/n_H
A       =10**4              #n_H/n_m

#define basic star stats, later will feed list ----->  UPDATE   THESE ARE TEMP VALUES NOW FOR TESTING DONT USE

#T        =5730.9728         #star temp in kelvin
#LogPe    =1.2594            #log_10 electron pressure in dynes cm-2
#theta    =5040.0/T

wavelength=5000.0           #angstroms

                            #calcultes ratio of ionized to neutroal hydrogen
def logsaha_h(T,LogPe):            #ie n+/n0, THESE SAY LOG BUT WILL ACTUALLY UNLOG THEM
    theta    =5040.0/T    
    logratio=-LogPe-13.595*theta+2.5*np.log10(T)-0.4772
    ratio=10**(logratio)
    return ratio

                            #calculates ratio of ionized to neutral metal
def logsaha_m(T,LogPe):            #ie m+/m0
    theta    =5040.0/T      
    logratio=-LogPe-7.9*theta+2.5*np.log10(T)-0.0971
    ratio=10**(logratio)    
    return ratio
    
def X(T,LogPe):                    #fraction H ionized
    ratio=logsaha_h(T,LogPe)
    x=ratio/(1+ratio)
    return x
                            #fraction metal ionized
def Y(T,LogPe):
    ratio=logsaha_m(T,LogPe)
    y=ratio/(1+ratio)
    return y
    
def gff(wavelength_angstroms,th):
    x=(1/wavelength_angstroms)*(1/1E-4)   #x is inverse wavelength in microns
    g=1.084+0.0188/th+(0.00161+0.02661/th)/x-(0.0192-0.03889/th+0.02833/th**2-0.007828/th**3+0.0007304/th**4)/x**2
    return g

def abc_m_index(m):  #this might have been the dumbest solution to this, oh well
    if m == 1:
        return [0.9916,0.09068,-0.2524]
    elif m == 2:
        return [1.105,-0.7922,0.4536]
    elif m == 3:
        return [1.101,-0.329,0.1152]
    elif m == 4:
        return [0.9736,0.0,0.0]
    elif m == 5:
        return [1.03,0.0,0.0]
    elif m == 6:
        return [1.097,0.0,0.0]
    elif m == 7:
        return [1.098,0.0,0.0]
    elif m > 7:
        return [1.0,0.0,0.0]
        

def gfb(wavelength_angstroms,m):
    x=(1.0/wavelength_angstroms)*(1.0/1.0E-4)   #x is inverse wavelength in microns
    a_m,b_m,c_m=abc_m_index(m)
    g=a_m+b_m/x+c_m/x**2
    return g
    
    
def u_m(T,m):
    u=(chi/(boltz*T))/m**2
    return u
    
    
#print X(),Y()
    
def alpha_lambda_H(T,LogPe):
    theta    =5040.0/T
    nu=speed_light/(wavelength*1.0E-8)
    x=(1.0/wavelength)*(1.0/1E-4)
    pt1= (2.0898E-14*e**(-u_m(T,1)))*(1-e**(-planck*nu/(boltz*T)))/((x**3.0)*2.0)
    pt2=0.0
    m0=1.0
    while u_m(T,m0) > planck*nu/(boltz*T):
        m0=m0+1.0 
    
    for mm in np.arange(m0,11,1):
        pt2=pt2+gfb(wavelength,mm)*e**(u_m(T,mm))/(mm**3.0)
    pt3=(1/(2*u_m(T,1)))*(e**(u_m(T,10))-1+gff(wavelength,theta))
    alpha=pt1*(pt2+pt3)
    return alpha

def kappa_lambda_H(T,LogPe):
    k=alpha_lambda_H(T,LogPe)*(1-X(T,LogPe))/((1+4*B)*m_H)
    return k
    
def alpha_lambda_H_ion(T,LogPe):  #man does this section look like a mess
    theta    =5040.0/T
    bigL=wavelength/1000.0
    nu=speed_light/(wavelength*1.0E-8)
    k_star=0.00680133 + 0.178708*bigL+0.164790*bigL**2 - 0.024842*bigL**3 + 5.95244E-4*bigL**4
    alphabf = 10**-26*(10**LogPe)*0.4158*theta**(5/2)*e**(1.726*theta)*(1-e**(-planck*nu/(boltz*T)))*k_star    
    alphaff = 10**-26 * 10**LogPe * (0.0053666 - 0.011493*theta + 0.027029*theta**2-(3.2062-11.924*theta+5.939*theta**2)*(wavelength/10**6)-(0.40192-7.0355*theta+0.34592*theta**2)*(wavelength**2/10**9))     
    alpha = alphabf+alphaff
    #print alphaff
    return alpha    
    
    
def kappa_lambda_H_ion(T,LogPe):
    k=alpha_lambda_H_ion(T,LogPe)*((1-X(T,LogPe))/((1+4*B)*m_H))
    return k

def sigma_r():
    r=5.799E-13/wavelength**4+1.422E-6/wavelength**6+2.784/wavelength**8
    return r

def sigma_R(T,LogPe):
    R=sigma_r()*(1-X(T,LogPe))/((1+4*B)*m_H)
    return R

def sigma_T(Temp,LogPe):
    T=sigma_t*(X(Temp,LogPe)+Y(Temp,LogPe)/A)/((1+4*B)*m_H) 
    return T

def sigma_lambda(T,LogPe):
    sig=sigma_R(T,LogPe)+sigma_T(T,LogPe) # you add the simgas together
    return sig    
    

    
    
def main():
    #k=kappa_lambda_H_ion()+kappa_lambda_H()
    data_out=[]
    for i in np.arange(0,len(data[:,0]),1):
        T=data[i,0]
        LogPe=data[i,1]
        LogH=np.log10(kappa_lambda_H(T,LogPe))
        LogHion=np.log10(kappa_lambda_H_ion(T,LogPe))
        LogSigma=np.log10(sigma_lambda(T,LogPe))
        data_out.append([T,LogPe,LogH,LogHion,LogSigma])          #CREATES A LINE WITH ALL THE CALCULATED VALUES AND ATTACHES IT TO THE END OF OUR DATA
    np.savetxt("data_out.dat",data_out,fmt='%1.3f')
    #print data_out
    diff= data_out-correct  # THIS CALCULATES THE DIFFERENCES BETWEEN MY DATA OUTPUT AND THE DATA OUTPUT ANA GAVE US, tells us how wrong we might be
    np.savetxt("diffs.dat",diff,fmt='%1.3f')
    #plotallthisshit(data_out)
    return data_out
    #print "Sum of kappas is "+ str(k)
    #print "Kappas before logs:"+str(kappa_lambda_H_ion())+" and "+str(kappa_lambda_H())
    #print "after logs:"+ str(np.log10(kappa_lambda_H_ion()))+" and " + str(np.log10(kappa_lambda_H()))



    
#kappa_lambda_H_ion()
    
# THIS SECTIONS RUNS THE CODE AND SPITS OUT SOME PLOTS, DIS IS DA END YO    
main()
data_out=np.array(main())   # UNFORTUNATELY YOU GOTTA TURN THIS THING BACK INTO A NUMPY ARRAY TO MAKE YOUR LIFE EASY
total=np.log10(10**data_out[:,2]+10**data_out[:,3]+10**data_out[:,4])
plt.plot(data_out[:,0],data_out[:,2],label="H",c='r',ls='-')
plt.plot(data_out[:,0],data_out[:,3],label="H-",c='b',ls='-')
plt.plot(data_out[:,0],data_out[:,4],label="sigma",c='g',ls='-')
plt.plot(data_out[:,0],total,label="total",c="cyan",ls='-')
plt.legend(bbox_to_anchor=(1, .41))
plt.ylabel("log opacity values (Kappa/gram)")       # LABELS ARE IMPORTANT
plt.xlabel("Temperature (K)")
plt.ion()
plt.show()   
