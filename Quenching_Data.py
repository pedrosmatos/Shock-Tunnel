import numpy as np

k = np.arange(18)
s = 5.0/3.0
fat = [s]
for i in np.arange(1,18): fat = np.append(fat, fat[i-1]*(s+i))
    

def sN2(T):
    c = [0.0, 0.88, 4.9, 48.0, 32.0]
    return  c[0] + c[1]*np.exp(-c[2]*300.0/T) + c[3]*np.exp(-c[4]*300.0/T)

def sO2(T):
    c = [25.1]
    return c[0]

def sNO(T):
    c = [43.0]
    return c[0]

def sO(T):
    c = [32.0]
    return c[0]

def sNO2(T):
    c = [82.0, 9.0, 0.54]
    return c[0] + c[1]*(300.0/T)**c[2]

def sNOBRE(T):
    c = [0.0]
    return c[0]

def sH2(T):
    c = [0.0]
    return c[0]

def sCO(T):
    c = [5.9,  5.3,  7.0,  22.1,  14.0] 
    return c[0] + c[1]*np.exp(-c[2]*300.0/T) + c[3]*np.exp(-c[4]*300.0/T)

def sH2O(T):
    c = [28.2,  3.39,  0.15,  2.95]
    x = c[2]*(300.0/T) + c[3]*(300.0/T)**2.0    
    S = ((x**s)*np.exp(-x)*(x**k)/fat).sum()    
    return c[0]*( (1.0+x)*np.exp(-x) + c[1]*S*x**(0.333) )

def sOH(T):
    c = [82.0]
    return c[0]

def sH(T):
    c = [12.0]
    return c[0]

def sN2O(T):
    c = [59.0,  0.99,  3.98,  0.16]
    x = c[2]*(300.0/T) +  c[3]*(300.0/T)**2.0
    S = ((x**s)*np.exp(-x)*(x**k)/fat).sum()
    return c[0]*( (1.0+x)*np.exp(-x) + c[1]*S*x**(0.333) )

def sCO2(T):
    c = [54.2,0.95,3.24,0.18]
    x = c[2]*(300.0/T) + c[3]*(300.0/T)**2.0
    S = ((x**s)*np.exp(-x)*(x**k)/fat).sum()
    return c[0]*( (1.0+x)*np.exp(-x) + c[1]*S*x**(0.333) )

def Fluorescence_Lifetime(P,T,X):
#X = xN2, xH2O, xCO2, xCO, xO2, xOH, xH, xO, xH2, xNO, xNOBRE, xN2O, xNO2
    X  = X/X.sum()
    kb = 1.38065E-23
    A  = 6.022141E23
    Ae = 4545454.5454

    sq = np.array([sN2(T),  sH2O(T),   sCO2(T),   sCO(T),
                   sO2(T),   sOH(T),     sH(T),    sO(T),
                   sH2(T),   sNO(T), sNOBRE(T),  sN2O(T),  
                   sNO2(T)])

    Mass = np.array([28.0,     18.0,      44.0,     28.0,
                     32.0,     17.0,       1.0,     16.0,
                      2.0,     30.0,      34.0,     44.0,
                     46.0])

    vNO = np.sqrt( 8.0*kb*T*A/(3.14159*(Mass[9]/1e3) ))

    Sq = (X*sq*(1.0 + Mass[9]/Mass)**(0.5)).sum()   
    Qm = 10**(-20.0)*vNO*(Sq)*P/(kb*T)
    return 10.0**9/(Qm+Ae)  