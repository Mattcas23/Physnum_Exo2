import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os

# Parameters
# TODO adapt to what you need (folder path executable input filename)
executable = r"/Users/a-x-3/Desktop/Ex2_2025_student/exe"  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"/Users/a-x-3/Desktop/Ex2_2025_student"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Name of the input file


nsteps = np.array([10, 20 , 50 , 100 , 150 , 200]) # TODO change
#nsteps = np.array([10])
nsimul = len(nsteps)  # Number of simulations to perform

tfin = 7776000  # Done : Verify that the value of tfin is EXACTLY the same as in the input file

dt = tfin / nsteps

energy = np.zeros(nsimul) # added 

paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

# Simulations
outputs = []  # List to store output file names
convergence_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    cmd = f"{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

#values = np.loadtxt("configuration.in.example")
delta = np.zeros(nsimul)

# Valeurs du fichier configutation.in.example. Vérifier à chaque fois si similaires

Omega = 0 
m = 0.075
L = 0.08
mu = 0.2
B0 = 0.01
w0 = np.sqrt( 12 * mu * B0 / ( m * L ** 2 ) )
print(w0)

theta0 = 1e-6
thetadot0 = 0.0

# ------------------------------ Valeurs Finales Pour Simulations -------------------------------#

for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation
    t = data[:, 0]

    theta = data[-1, 1]  # final position, velocity, energy
    thetadot = data[-1, 2]
    emec = data[-1, 3]
    pnc = data[-1, 4]
    convergence_list.append(theta)
    
    theta_an = theta0 * np.cos(w0 * data[-1,0])  # solution analytique pour la position 
    thetadot_an = - theta0 * np.sin(w0 * data[-1,0]) * w0 # solution analytique pour la vitesse
    print (f"theta_an : {theta_an}")
    print (f"thetadot_an : {thetadot_an}")
    print (f"theta : {theta}")
    print (f"thetadot : {thetadot}")

    delta[i] =  np.sqrt( w0**2 * (theta - theta_an)**2 + (thetadot - thetadot_an)**2 ) # delta simulation 
    # TODO compute the error for each simulation

lw = 1.5
fs = 16

# ------------------------------------------- Zone Figures -------------------------------------- # 

def Emec () : # Affiche l'Emec en fonction du temps

    fig, ax = plt.subplots(constrained_layout=True)
    plt.ticklabel_format(axis='y', style='scientific', scilimits = (-3,-3))
    ax.plot(t, data[:,3], color = 'orange' , label = '$n_{step} = $' + f"{nsteps[-1]:.0f}")
    ax.plot(t, np.mean(data[:,3])*np.ones(t.size) , color = "grey" , linestyle = "dashed" , label = '$<E_{mec}> = $' + f"{np.mean(data[:,3]):.4e}" ) 
    ax.set_ylabel('$E_{mec}$', fontsize=fs)
    ax.set_xlabel('$\\Delta t$ [s]', fontsize=fs)
    plt.legend()

def Pnc () : # Affiche la puissance des forces non-conservatives en fonction du temps (ajouter la dérivée de Emec après) 

    fig, ax = plt.subplots(constrained_layout=True)
    plt.ticklabel_format(axis='y', style='scientific', scilimits = (-3,-3))
    ax.plot(t, data[:,4], color = 'orange' , label = '$n_{step} = $' + f"{nsteps[-1]:.0f}")
    ax.set_ylabel('$P_{nc}$', fontsize=fs)
    ax.set_xlabel('$\\Delta t$ [s]', fontsize=fs)
    plt.legend() 

def Theta () : # Affiche la position en fonction du temps 
    
    fig, ax = plt.subplots(constrained_layout=True)
    ax.plot(t, 1e-06 * np.cos(w0 * t) , color = 'red' , label = '$n_{step} = $' + f"{nsteps[-1]:.0f}")    
    ax.plot(t, data[:,1], color = 'black' , label = '$n_{step} = $' + f"{nsteps[-1]:.0f}", linestyle = "dashed")
    ax.set_ylabel('$\\theta$', fontsize=fs)
    ax.set_xlabel('$\\Delta t$ [s]', fontsize=fs)
    plt.legend()

def Thetadot () : # Affiche la vitesse en fonction du temps 

    fig, ax = plt.subplots(constrained_layout=True)
    #ax.plot(t, data[:,1], color = 'blue' , label = '$n_{step} = $' + f"{nsteps[-1]:.0f}")
    ax.plot(t, - theta0 * np.sin(w0 * t) * w0, color = 'red' , label = 'Analytique : $n_{step} = $' + f"{nsteps[-1]:.0f}")
    ax.plot(t, data[:,2], color = 'black' , label = '$n_{step} = $' + f"{nsteps[-1]:.0f}", linestyle = "dashed")
    ax.set_ylabel('$\\dot{\\theta}$', fontsize=fs)
    ax.set_xlabel('$\\Delta t$ [s]', fontsize=fs)
    plt.legend()


def Delta (n_order = 1) : # Affiche l'erreur en fonction du nombre de pas 

    plt.figure()
    plt.loglog(1/nsteps[:-1], delta[:-1], 'r+-', linewidth=lw)
    #plt.loglog(nsteps, pow(nsteps,n_order), color = 'black' ,linewidth = lw , label = f"$1/N^{n_order}$" , linestyle = 'dashed')
    plt.xlabel('$\\Delta t$', fontsize=fs)
    plt.ylabel('$\\delta$', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.grid(True)
    plt.legend()

def Phase () :

    fig, ax = plt.subplots(constrained_layout=True)
    ax.plot(data[:,1], data[:,2], color = 'red' , label = '$n_{step} = $' + f"{nsteps[-1]:.0f}")
    ax.set_ylabel('$\\dot{\\theta}$', fontsize=fs)
    ax.set_xlabel('$\\theta$', fontsize=fs)
    plt.legend()

def PointCarre () :


    print("A faire")


# ------------------ Affichage ------------------ # 

#Emec ()
#Pnc ()
Delta ()
Theta () 
Thetadot()
Phase()
#PointCarre () 

plt.show()
