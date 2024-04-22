import matplotlib as mpl

fsize = 18
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams["figure.figsize"] = (6,5)
mpl.rcParams['axes.labelsize'] = fsize
mpl.rcParams['xtick.labelsize'] = fsize
mpl.rcParams['ytick.labelsize'] = fsize
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['mathtext.fontset'] = 'dejavusans'
mpl.rcParams.update({'font.size': fsize})

import matplotlib.pyplot as plt

import sympy as sp
import numpy as np
import scipy.integrate as integrate

#Coulomb barrier

Z_9Be = 4
Z_n = 1

A_9Be = 9
A_n = 1

R_9Be = 1.2*A_9Be**(1/3) #fm
R_n = 0.8 #fm

e2_4pi_e0 = 1.44 #Mev/fm

E_CB = Z_9Be*Z_n*e2_4pi_e0/(R_9Be+R_n)

print(E_CB)

m_9Be = 9.01218306 #u
m_n = 1.00866491595 #u

#1u = 931.5 MeV/c^2
mu_c2 = (m_9Be*m_n)*931.5/(m_9Be+m_n) #MeV
print(mu_c2)

hb_c = 197.3269804 #MeV fm

fig_e, ax_e = plt.subplots()
fig_a, ax_a = plt.subplots(figsize=(6, 4.9))
fig_v, ax_v = plt.subplots()

def sin_d_l(k,a):
    return (2*mu_c2/hb_c**2)*np.log((4*k**2)/(a**2) + 1)/(4*k)

# Define the function
def integrand(x, k, a):
    return  (np.sin(k*x))**2*np.exp(-a*x)/(k*x)

E_values = np.arange(np.ceil(E_CB), 101, 1) #MeV
E_values_V = np.arange(10, 101, 10) #MeV

a_values = np.arange(1, 50, 1)
a_values_V = np.arange(5, 51, 5)
V0_values = np.arange(1, 101, 1) #MeV

dl_values_E = [[0]*len(E_values)]*len(a_values_V)
dl_values_a = [[0]*len(a_values)]*len(E_values_V)
dl_values_v = [[0]*len(V0_values)]*len(E_values_V)

V0_val = 50 #MeV
for a, a_val in enumerate(a_values_V):
    for e, e_val in enumerate(E_values):
        k_val = np.sqrt(2*mu_c2*e_val)/hb_c

        dl_values_E[a][e] = np.degrees(np.arcsin(V0_val*sin_d_l(k_val, a_val)))

    ax_e.plot(E_values, dl_values_E[a], label=str(a_val))

for e, e_val in enumerate(E_values_V):
    k_val = np.sqrt(2*mu_c2*e_val)/hb_c
    for a, a_val in enumerate(a_values):
        dl_values_a[e][a] = np.degrees(np.arcsin(V0_val*sin_d_l(k_val, a_val)))
    
    ax_a.plot(a_values, dl_values_a[e], label=str(e_val))

a_val = 25 #1/fm
for e, e_val in enumerate(E_values_V):
    k_val = np.sqrt(2*mu_c2*e_val)/hb_c
    for v0, V0_val in enumerate(V0_values):
        dl_values_v[e][v0] = np.degrees(np.arcsin(V0_val*sin_d_l(k_val, a_val)))
    
    ax_v.plot(V0_values, dl_values_v[e], label=str(e_val))

ax_e.grid()
#ax_e.set_title(label="Energy deposition")
ax_e.set_xlabel("Energy [MeV]")
ax_e.set_ylabel(r"$\delta_0$ [$^\circ$]")
ax_e.set_yscale(value="log")

ax_a.grid()
#ax_a.set_title("Photon count")
ax_a.set_xlabel(r"$\alpha$ [fm$^{-1}$]")
ax_a.set_ylabel(r"$\delta_0$ [$^\circ$]")
ax_a.set_yscale(value="log")

ax_v.grid()
#ax_v.set_title("Photon count")
ax_v.set_xlabel(r"V$_0$ [MeV fm]")
ax_v.set_ylabel(r"$\delta_0$ [$^\circ$]")

ax_e.text(20, 0.015, r"$V_0=$ 50 MeV")
ax_a.text(10, 40, r"$V_0=$ 50 MeV")
ax_v.text(30, 0.8, r"$\alpha=25$ fm$^{-1}$")

ax_e.legend(loc="lower right", title="a [1/fm]", title_fontsize=12)
ax_a.legend(title="E [MeV]", title_fontsize=12)
ax_v.legend(title="E [MeV]", title_fontsize=12)

fig_e.savefig("phase_shift_vs_E.pdf", bbox_inches="tight")
fig_a.savefig("phase_shift_vs_alpha.pdf", bbox_inches="tight")
fig_v.savefig("phase_shift_vs_V0.pdf", bbox_inches="tight")