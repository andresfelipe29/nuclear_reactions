import matplotlib.pyplot as plt
import numpy as np

mt = 63.9291418 #[AMU]
mp = 9.01212306 #[AMU]

Elab = 19 #[MeV] #[AMU]

M = mt+mp #[AMU]
miu = (mt*mp)/(mt+mp) #[AMU]

Ecm = (mt/M)*Elab

print(Ecm)

Lambda = mp/mt

print("Lambda=", Lambda)

DO_n_J = 2.9e27 #[cm-2]

e2_4piE0 = 14.4e-14 # [MeV*cm]

Zp = 4
Zt = 30

a = 0.5*(e2_4piE0*Zp*Zt/Ecm)

alpha = np.radians([30, 60, 90, 120, 150])
N_al = np.array([32983, 2296, 545, 181, 82])
N_al_err = np.array([181, 50, 10, 7, 3])

def cos_t(a):
    return -Lambda*(np.sin(a)**2) + (np.cos(a)*np.sqrt(1-(Lambda**2)*(np.sin(a)**2)))

def calc_cross_sec_lab(n): # [b]
    cross_sec = n/DO_n_J*1e24
    cross_sec_err = N_al_err/DO_n_J*1e24
    return cross_sec, cross_sec_err 

def cal_cross_sec_cm(cross_sec_lab, corss_sec_lab_err, cos_t): # [b]
    p = ((np.abs(1+Lambda*cos_t))/np.power(1+2*Lambda*cos_t+Lambda**2, 1.5))
    cross_sec = cross_sec_lab*p
    cross_sec_err = corss_sec_lab_err*p

    return cross_sec, cross_sec_err

def Rutherford(t): # [b]
    return (a**2/4)/(np.sin(t/2)**4)*1e24

cos_t_values = cos_t(alpha)
print("cos_t_values")
print(cos_t_values)

print("t_values")
print(np.degrees(np.arccos(cos_t_values)))

cross_sec_lab, corss_sec_lab_err = calc_cross_sec_lab(N_al)
print("lab")
print(cross_sec_lab, corss_sec_lab_err)

cross_sec_cm, cross_sec_cm_err = cal_cross_sec_cm(cross_sec_lab, corss_sec_lab_err, cos_t_values)
print("cm")
print(cross_sec_cm, cross_sec_cm_err)

x = np.radians(np.arange(30, 161, 1))

plt.grid()

plt.errorbar(x=np.degrees(np.arccos(cos_t_values)), y=cross_sec_cm, yerr=cross_sec_cm_err, label="C.M", fmt="o", ms=3)
plt.plot(np.degrees(x), Rutherford(x), label="Ruth")

plt.xlabel(r"$\theta$ [$^\circ$]")
plt.ylabel(r"$d\sigma/d\Omega$ [b]")

plt.yscale(value="log")

plt.legend()
plt.savefig("seccion_eficaz_cm.pdf", bbox_inches="tight")

plt.clf()

x = np.arccos(cos_t_values)

ruth_t_values = Rutherford(x)

ratio = cross_sec_cm/ruth_t_values

ratio_err = cross_sec_cm_err/ruth_t_values

plt.grid()

plt.errorbar(x=np.degrees(x), y=ratio, yerr=ratio_err, label="C.M/Ruth", fmt="o", ms=2)

plt.xlabel(r"$\theta$ [$^\circ$]")
plt.ylabel("C.M/Rutherford")

plt.legend()

plt.plot()
plt.savefig("seccion_eficaz_ratio.pdf", bbox_inches="tight")
