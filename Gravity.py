import numpy as np
# for graphing
import matplotlib.pyplot as plt
from  matplotlib import animation
from IPython.display import HTML
plt.rcParams['figure.figsize']=[22,22]
plt.rcParams.update({'font.size':16})
plt.rcParams['animation.ffmpeg_path'] = r'C:\Users\jaket\Desktop\ffmpeg-2021-06-27-git-49e3a8165c-full_build\bin\ffmpeg.exe'
#You will need to change above to your ffmpeg path (may be able to comment out above line on certain broswers/OS)
plt.rcParams['animation.embed_limit'] = 20971520.0*200 #about 400 MB
Writer = animation.writers['ffmpeg']
writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)
import timeit

class Gravity:

    def __init__(self, rs: dict, vs: dict, ms: list) -> None:
        self.G = 6.67*10**(-11)
        #define initial conditions
        self.rs = rs
        self.vs = vs
        self.ms = ms
        self.run = False #tell whether the system has been simulated some amount
        self.n = len(ms) #number of bodies




        


    def force(self,x1,x2, m1: float, m2: float, G: float): #force on block x1
        denom = (np.linalg.norm(x2-x1)**3)
        if denom == 0:
            denom += 10**(-6)

        prefactor = m1*m2*G/denom
        return prefactor*(x2-x1)


    def calculate_Trajectory(self, t: float, iters: int): #input dictionary of initial coordinates and velocities
        interval = t/iters #time resolution
        rvals = {}
        vvals = {}
        forces = {}
        #let's generate the dictionaries of all values at all timesteps

        self.interval = interval
        self.t = t
        self.iters = iters
        

        #first index = timestep ======= in this program index i will be denoting the timestep
        #second index = what body is this a property of/on = index j denotes which body

        #initialize all of the forces to zero
        for i in range(0, iters+1, 1):
            forces[i] = {}
            vvals[i] = {}
            rvals[i] = {}
            for j in range(self.n):
                forces[i][j] = np.array([0,0,0], dtype = 'float64')
                rvals[i][j] = np.array([0,0,0], dtype = 'float64')
                vvals[i][j] = np.array([0,0,0], dtype = 'float64')

        #initialize the zeroth timestep
        rvals[0] = self.rs
        vvals[0] = self.vs
        for i in range(self.n): #here i is the first body, j is the second body in the pairwise interaction
            rvals[0][i] = np.array(rvals[0][i], dtype = 'float64')
            vvals[0][i] = np.array(vvals[0][i], dtype = 'float64')
            for j in range(i+1,self.n,1):
                print(i, j)
                forces[0][i] += self.force(x1 = rvals[0][i], x2 = rvals[0][j], m1 = self.ms[i], m2 = self.ms[j], G = self.G)
                forces[0][j] -= forces[0][i]
                
        #now loop through and generate all of the time steps
        for i in range(1, iters+1, 1):
            for j in range(self.n):
                for k in range(j+1, self.n, 1):
                    forces[i][j] += self.force(rvals[i-1][j], rvals[i-1][k], self.ms[j], self.ms[k], self.G)
                    forces[i][k] -= self.force(rvals[i-1][j], rvals[i-1][k], self.ms[j], self.ms[k], self.G)
            
            for j in range(self.n):   
                vvals[i][j] = vvals[i-1][j] + forces[i][j]*interval/(ms[j]) #velo at the ith timestep on the jth body
                rvals[i][j] = rvals[i-1][j] + vvals[i][j]*interval
            #print(i*interval, 'Timestep Done')
        self.rout = rvals
        self.vout = vvals
        self.fout = forces

    def display_Trajectory(self):
        fig = plt.figure(figsize=(22,22))
        ax = plt.axes(projection='3d')
        ax.view_init(elev = 90, azim = 0)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        colors = plt.cm.rainbow(np.linspace(0,1,self.n))


        for j in range(self.n):
            xvals = []
            yvals = []
            zvals = []
            
            
            for i in range(len(self.rout)):
                xvals.append(self.rout[i][j][0])
                yvals.append(self.rout[i][j][1])
                zvals.append(self.rout[i][j][2])
            ax.plot(xvals, yvals, zvals, color = colors[j], ls = '--', marker = '.', label = f'Body {j+1}')
        ax.legend()
        plt.show()

    def animate_it(self):
        fig = plt.figure(figsize=(22,22))
        ax = plt.axes(projection='3d')
        ax.view_init(azim=45, elev=6)

        
        ax.set_xlim(-2,2)
        ax.set_ylim(-2,2)
        ax.set_zlim(-2,2)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        colors = plt.cm.rainbow(np.linspace(0,1,self.n)) 

        #initialize the scatter plot
        plot = [ax.scatter(self.rout[0][i][0], self.rout[0][i][1], self.rout[0][i][2], marker = 'o', color = colors[i], label = f'Body {i+1}') for i in range(self.n)]

        ax.legend()

        def animate(i, rvals, plot):
            ax.cla()
            ax.set_xlim(-2,2)
            ax.set_ylim(-2,2)
            ax.set_zlim(-2,2)
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plot = [ax.scatter(rvals[i][j][0], rvals[i][j][1], rvals[i][j][2], marker = 'o', color = colors[j], label = f'Body {j}') for j in range(self.n)]
            ax.legend()
            return plot
        


        ani = animation.FuncAnimation(fig = fig, func = animate, frames = self.iters+1, fargs=(self.rout, plot),
                                        interval=self.interval/8, blit=False, repeat=False)

        ani.save('3bodysim.mp4', writer = writer)
        #HTML(ani.to_jshtml())
        #plt.show()

       



if __name__ == '__main__':
    
    
    G = 6.67*10**(-11) #this is used to calculate ideal times and speeds for clean orbits

    ms = [10**15, 1, 1, 1]
    velo_test = np.sqrt(G*ms[0]/2) #should yield a clean orbit for a while


    rs = {
        0: [0,0,0],
        1: [0,2,0],
        2: [2,0,0],
        3: [-2,0,0]
        }
    vs = {
        0: [0, 0, 0],
        1: [0, 0, velo_test],
        2: [0, 0, velo_test],
        3: [0, velo_test, 0]
    }

    t = 2*np.pi*2*np.sqrt(2/(G*10**(15)))*1.25
    iters = 1200
    a = Gravity(rs, vs, ms)
    a.calculate_Trajectory(t, iters)
    #a.display_Trajectory()
    a.animate_it()
