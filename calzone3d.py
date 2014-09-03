from numpy import *
from scipy.spatial import Delaunay
from scipy.special import cbrt
from itertools import combinations
from math import pi
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from collections import defaultdict
from matplotlib.patches import Rectangle, Ellipse



#### Calzone Constants
K = {'sc':0.01, 'afzy':1.0, 'sxm':0.0005, 'astgp':0.0, 'in':0.15,
     'dcp':0.01, 'ifzy':0.2, 'sxp':0.001, 'astgpp':1.0, 'out':0.04,
     'dnp':0.01, 'weep':0.005, 'dmp':0.002, 'istg':0.3, 'ins':0.08,
     'dnpp':1.5, 'weepp':1.0, 'dmpp':0.2, 'awee':0.3, 'outs':0.02,
     'aie':1.0, 'stgp':0.2, 'sstg':0.02, 'iweep':0.01, 'inw':0.04,
     'iie':0.4, 'stgpp':2.0, 'dstg':0.015, 'iweepp':1.0, 'outw':0.01,
     'ez':0.5}

koutarr = array([K['out'], K['out'], K['outs'], K['outs'], K['outw'], K['outw']])
kinarr = array([K['in'], K['in'], K['ins'], K['ins'], K['inw'], K['inw']])

J = {'istg':0.05, 'aie':0.01, 'awee':0.05, 'iie':0.01, 'iwee':0.05,
     'afzy':0.01, 'm':0.05, 'ifzy':0.01, 'astg':0.05}

eps = 0.00007
N = 1.0
Wee1T = 0.8

namesN = {'MPF':0, 'preMPF':1, 'IE':2, 'FZY':3, 'StgP':4, 'Stg':5, 'Wee1P':6, 'Wee1':7}
namesC = {'MPF':0, 'preMPF':1, 'StgP':2, 'Stg':3, 'Wee1P':4, 'Wee1':5}
transport_array = array([0, 1, 4, 5, 6, 7])



#### Physical / Geometric Constants
gw, gh, gd = (20, 10, 10)
a, b, c = (18, 8, 8)
aa, bb, cc = (float(a)*a, float(b)*b, float(c)*c)
volumeT = 4.0*pi*a*b*c / 3.0
volumeN = eps * volumeT
radiusN = cbrt((3*volumeN) / (4*pi))

D = 0.01
boundary_force = 250.0
nuclear_force = 50.0
nuclear_distance = 1.0
visc = 0.85

h = min(0.2/D, 0.01)
hh = h*h
r = D*h
rr = 1.0 - 6*D*h



#### Create Ellipse Mask
grid = zeros((2*gw + 1, 2*gh + 1, 2*gd + 1))
for i in range(2*gw + 1):
    for j in range(2*gh + 1):
        for k in range(2*gd + 1):
            x, y, z = (i - gw, j - gh, k - gd)
            grid[i,j,k] = x*x/aa + y*y/bb + z*z/cc
interior_mask = grid < 1.0
exterior_mask = ~interior_mask
embryo_depth = interior_mask.sum(2) + 1.0


#### Initialize Structures
Y = zeros((1, 8))
X = zeros((1, 3))
V = zeros((1, 3))
A = zeros((1, 3))

##W = ones((2*gw + 1, 2*gh + 1, 2*gd + 1, 6))
##W = random.rand(2*gw + 1, 2*gh + 1, 2*gd + 1, 6)
W = zeros((2*gw + 1, 2*gh + 1, 2*gd + 1, 6))

initial_values = array([[1.0, 0, 0.8, 0, 0.8, 0]])
W[:,:,:] = initial_values
W[exterior_mask] = zeros(W[exterior_mask].shape)
mask_shape = W[exterior_mask].shape




def draw_embryo(im_num, spec1 = None, spec2 = None):
    
    
    ax = plt.gca()

    plt.cla()
    if spec1:
        alphas = W.sum(2)[:,:,namesC[spec1]]
        trans = (0.5 if spec2 else 1.0)
        plt.imshow(alphas.transpose()/embryo_depth.transpose(), cmap = cm.Oranges, vmin = 0, vmax = 1, alpha = trans)
    if spec2:
        alphas = W.sum(2)[:,:,namesC[spec2]]
        plt.imshow(alphas.transpose()/embryo_depth.transpose(), cmap = cm.PuBu, vmin = 0, vmax = 1, alpha = 0.5)

    for n in range(N):
        ax.add_patch(Ellipse((X[n,0] + gw, X[n,1] + gh), radiusN, radiusN))
##    plt.plot(X[:,0] + gw, Y[:,1] + gh, 'bo')

    plt.text(0.5, 0.5, str(t))

    plt.savefig('im/' + str(im_num) + '.png', bbox_inches = 'tight')


#### Simulation Loop
t = 0.0
stop_t = 1.0
steps = int(stop_t/h)
T = []

record1 = []
record2 = []
record3 = []
record4 = []

print 'Begin Simulation'
print 'Total time:', stop_t
print 'Step size:', h
print 'Steps:', steps

for step in range(steps):

    N = X.shape[0]

    #########################################################################################
    #### Nuclei Movement ####################################################################

    #### End Nuclei Movement ################################################################
    #########################################################################################


    

    #########################################################################################
    #### Diffusion ##########################################################################
    
    new_W = zeros((2*gw + 1, 2*gh + 1, 2*gd + 1, 6))

    # look forward
    new_W[1:-1,1:-1,1:-1] += r * W[0:-2,1:-1,1:-1]
    # reflect back
    new_W[1:-1,1:-1,1:-1] += new_W[2:,1:-1,1:-1] * exterior_mask[2:,1:-1,1:-1, newaxis]
    new_W[exterior_mask] = zeros(mask_shape)

    # look backward
    new_W[1:-1,1:-1,1:-1] += r * W[2:,1:-1,1:-1]
    # reflect back
    new_W[1:-1,1:-1,1:-1] += new_W[0:-2,1:-1,1:-1] * exterior_mask[0:-2,1:-1,1:-1, newaxis]
    new_W[exterior_mask] = zeros(mask_shape)

    # look up
    new_W[1:-1,1:-1,1:-1] += r * W[1:-1,0:-2,1:-1]
    # reflect back
    new_W[1:-1,1:-1,1:-1] += new_W[1:-1,2:,1:-1] * exterior_mask[1:-1,2:,1:-1, newaxis]
    new_W[exterior_mask] = zeros(mask_shape)

    # look down
    new_W[1:-1,1:-1,1:-1] += r * W[1:-1,2:,1:-1]
    # reflect back
    new_W[1:-1,1:-1,1:-1] += new_W[1:-1,0:-2,1:-1] * exterior_mask[1:-1,0:-2,1:-1, newaxis]
    new_W[exterior_mask] = zeros(mask_shape)

    # look left
    new_W[1:-1,1:-1,1:-1] += r * W[1:-1,1:-1,0:-2]
    # reflect back
    new_W[1:-1,1:-1,1:-1] += new_W[1:-1,1:-1,2:] * exterior_mask[1:-1,1:-1,2:, newaxis]
    new_W[exterior_mask] = zeros(mask_shape)

    # look right
    new_W[1:-1,1:-1,1:-1] += r * W[1:-1,1:-1,2:]
    # reflect back
    new_W[1:-1,1:-1,1:-1] += new_W[1:-1,1:-1,0:-2] * exterior_mask[1:-1,1:-1,0:-2, newaxis]
    new_W[exterior_mask] = zeros(mask_shape)

    W[1:-1,1:-1,1:-1] = rr * W[1:-1,1:-1,1:-1] + new_W[1:-1,1:-1,1:-1]
    
    #### End Diffusion ######################################################################
    #########################################################################################


    new_W = copy(W)
    new_Y = copy(Y)


    #########################################################################################
    #### Transport ##########################################################################

    intersections = defaultdict(dict)
    for i in range(N):
        x, y, z = X[i, :]
        coords = (int(floor(x)+gw), int(floor(y)+gh), int(floor(z)+gd))
        intersections[i][coords] = 1.0
        
    for n in intersections.keys():
        for voxel, fraction in intersections[n].iteritems():
            dilution = fraction*volumeN / (1.0 - fraction*volumeN)
            new_W[voxel] += fraction*koutarr*Y[n, transport_array]*dilution  ## n-->c
            new_W[voxel] -= fraction*kinarr*new_W[voxel]*dilution  ## c-->n
            new_Y[n][transport_array] -= fraction*koutarr*Y[n, transport_array]  ## n-->c
            new_Y[n][transport_array] += fraction*kinarr*W[voxel]   ## c-->n            
            
    #### End Transport ######################################################################
    #########################################################################################




    #########################################################################################
    #### Reaction - Cytoplasm ###############################################################

    new_W[:,:,:,0] += h * (K['sc'] - (K['dcp'] + K['weep'] + K['weepp']*W[:,:,:,5])*W[:,:,:,0] + (K['stgp'] + K['stgpp']*W[:,:,:,2])*W[:,:,:,1])
    new_W[:,:,:,1] += h * -((K['dcp'] + K['stgp'] + K['stgpp']*W[:,:,:,2])*W[:,:,:,1] + (K['weep'] + K['weepp']*W[:,:,:,5])*W[:,:,:,0])
    new_W[:,:,:,2] += h * -(K['dstg']*W[:,:,:,2] + ((K['astgp'] + K['astgpp']*W[:,:,:,0])*W[:,:,:,3])/(J['astg'] + W[:,:,:,3]) - (K['istg']*W[:,:,:,2])/(J['istg']+W[:,:,:,2]))
    new_W[:,:,:,3] += h * (K['sstg'] - K['dstg']*W[:,:,:,3] - ((K['astgp'] + K['astgpp']*W[:,:,:,0])*W[:,:,:,3])/(J['astg'] + W[:,:,:,3]) + (K['istg']*W[:,:,:,2])/(J['istg']+W[:,:,:,2]))
    new_W[:,:,:,4] += h * (-(K['awee']*W[:,:,:,4])/(J['awee'] + W[:,:,:,4]) + ((K['iweep'] + K['iweepp']*W[:,:,:,0])*W[:,:,:,5])/(J['iwee'] + W[:,:,:,5])) ## should wee1pc be determined using wee1t?
    new_W[:,:,:,5] += h * ((K['awee']*W[:,:,:,4])/(J['awee'] + W[:,:,:,4]) - ((K['iweep'] + K['iweepp']*W[:,:,:,0])*W[:,:,:,5])/(J['iwee'] + W[:,:,:,5]))

    #### End Reaction - Cytoplasm ###########################################################
    #########################################################################################
    



    #########################################################################################
    #### Reaction - Nuclei ##################################################################

    new_Y[:,0] += h * (-(K['dnp'] + K['dnpp']*Y[:,3] + K['weep'] + K['weepp']*Y[:,7])*Y[:,0] + (K['stgp'] + K['stgpp']*Y[:,4])*Y[:,1])
    new_Y[:,1] += h * (-(K['dnp'] + K['dnpp']*Y[:,3] + K['stgp'] + K['stgpp']*Y[:,4])*Y[:,1] + (K['weep'] + K['weepp']*Y[:,7])*Y[:,0])
    new_Y[:,2] += h * ((K['aie']*Y[:,0]*(1.0 - Y[:,2]))/(J['aie'] + 1.0 - Y[:,2]) - (K['iie']*Y[:,2])/(J['iie'] + Y[:,2]))
    new_Y[:,3] += h * ((K['afzy']*Y[:,2]*(1.0 - Y[:,3]))/(J['afzy'] + 1.0 - Y[:,3]) - (K['ifzy']*Y[:,3])/(J['ifzy'] + Y[:,3]))
    new_Y[:,4] += h * (-K['dstg']*Y[:,4] + ((K['astgp'] + K['astgpp']*Y[:,0])*Y[:,5])/(J['astg'] + Y[:,5]) - (K['istg']*Y[:,4])/(J['istg'] + Y[:,4]))
    new_Y[:,5] += h * (-K['dstg']*Y[:,5] - ((K['astgp'] + K['astgpp']*Y[:,0])*Y[:,5])/(J['astg'] + Y[:,5]) + (K['istg']*Y[:,4])/(J['istg'] + Y[:,4]))
    new_Y[:,6] += h * (-(K['awee']*Y[:,6])/(J['awee'] + Y[:,6]) + ((K['iweep'] + K['iweepp']*Y[:,0])*Y[:,7])/(J['iwee'] + Y[:,7]))
    new_Y[:,7] += h * ((K['awee']*Y[:,6])/(J['awee'] + Y[:,6]) - ((K['iweep'] + K['iweepp']*Y[:,0])*Y[:,7])/(J['iwee'] + Y[:,7]))

    #### End Reaction - Nuclei ##############################################################
    #########################################################################################

##    new_W[:,:,:,4] = Wee1T - new_W[:,:,:,5]

    W = copy(new_W)
    W[exterior_mask] = zeros(mask_shape)
    Y = copy(new_Y)

   

    if not step % 10:
        print step
        draw_embryo(step, 'MPF')

    record1.append(W[gw,gh,gd])
    record2.append(Y[0])
##    record3.append(Y[1])
    T.append(t)
    t += h
    if t > stop_t:
        break



if record1:
    plt.cla()
    if record2:
        if record3:
            plt.subplot(221)
            plt.plot(T, record1)
            plt.subplot(222)
            plt.plot(T, record2)
            plt.subplot(223)
            plt.plot(T, record3)
            if record4:
                plt.subplot(224)
                plt.plot(T, record4)
        else:
            plt.subplot(211)
            plt.plot(T, record1)
            plt.subplot(212)
            plt.plot(T, record2)
    else:
        plt.plot(T, record1)
    plt.show()



















