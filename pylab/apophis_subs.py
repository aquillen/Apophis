
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import matplotlib.tri as tri

from scipy.spatial.transform import Rotation as R  # for rotations
from scipy import linalg   # for eigenvalues 

hour = 60*60.


# read in springs file
def read_springs(mdir,froot,index):
    #springs[i].i, springs[i].j,
    # springs[i].ks, springs[i].rs0, springs[i].gamma, springs[i].k_heat
    sfname = mdir + froot + '_'+ iostring(index) + '_springs.txt'
    print(sfname)
    mtab = np.loadtxt(sfname)
    i_index = np.array(mtab[:,0], dtype=int)
    j_index = np.array(mtab[:,1], dtype=int)
    ks = mtab[0,2]  #all the same
    rs0 = mtab[:,3] # rest length, vector 
    ggamma = mtab[0,4] # all same 
    return i_index,j_index,ks,rs0,ggamma
    
# read in particles file
def read_particles(mdir,froot,index):
    # r->particles[i].x, r->particles[i].y, r->particles[i].z,
    #     r->particles[i].vx, r->particles[i].vy, r->particles[i].vz,
    #     r->particles[i].r, r->particles[i].m);
    sfname = mdir + froot + '_'+ iostring(index) + '_particles.txt'
    print(sfname)
    mtab = np.loadtxt(sfname,skiprows=1)
    xarr = mtab[:,0]
    yarr = mtab[:,1]
    zarr = mtab[:,2]
    vxarr = mtab[:,3]
    vyarr = mtab[:,4]
    vzarr = mtab[:,5]
    marr  = mtab[:,7]  
    return xarr,yarr,zarr,vxarr,vyarr,vzarr,marr

# file number manipulation
def iostring(index):
    junks = ''
    if (index < 10):
        junks += '0'
    if (index < 100):
        junks += '0'
    if (index < 1000):
        junks += '0'
    if (index < 10000):
        junks += '0'
    if (index < 100000):
        junks += '0'
    tstring = '{:d}'.format(index)
    junks += tstring
    #print(junks)
    return junks
    
# given a list of eigenvalues and eigenvectors of a 3x3 matrix, return a sorted list
# smallest eigenvalue first 
# also take real parts (assuming via a symmetric matrix)
def sort_eigs_3(w,v):
    jsort = np.argsort(w) # arguments/indices of a sorted array of eigenvalues, 
    # low to high
    jmax = jsort[2]  # index of maximum eigenvalue
    vmax = np.squeeze(np.asarray(v[:,jmax]))   # corresponding eigenvector
    jmin = jsort[0]
    vmin = np.squeeze(np.asarray(v[:,jmin]))   
    jmid = jsort[1]
    vmid = np.squeeze(np.asarray(v[:,jmid]))   
    eigmax = np.real(w[jmax])  #largest eigenvalue 
    eigmid = np.real(w[jmid])
    eigmin = np.real(w[jmin])
    wsort = np.array([eigmin,eigmid,eigmax])  # eigenvalues order smallest to largest  
    vsort = np.array([vmin,vmid,vmax])  # eigenvectors in same order 
    det_rot = np.linalg.det(vsort)
    # print('det_rot ',det_rot) # half the time we had -1 value 
    if (det_rot <0):
        vsort = np.array([-vmin,vmid,vmax])  # flip one so det becomes positive?
         
    return wsort, vsort 


# read in ext body output file to return spring potential and rotational rotational energy
# t x y z vx vy vz omx omy omz llx lly llz Ixx Iyy Izz Ixy Iyz Ixz KErot PEspr PEgrav Etot dEdtnow
# 0 1 2 3 4  5  6  7  8     9  10  11  12  13   14  15 16   17 18  19    20
# read ext
def read_ext(mdir,froot):
    #springs[i].i, springs[i].j,
    # springs[i].ks, springs[i].rs0, springs[i].gamma, springs[i].k_heat
    sfname = mdir + froot  + '_ext.txt'
    print(sfname)
    mtab = np.loadtxt(sfname)
    KErot = np.array(mtab[:,19])
    PEspr = np.array(mtab[:,20])
    tarr = np.array(mtab[:,0])
    return tarr,KErot,PEspr


# store all simulation output info in a single class structure 
class sim_struct():
    def __init__(self,mdir,froot,index,shape_vol):
        self.index = index  # which output
        self.mdir  = mdir  # directory 
        self.froot = froot  # fileroot!
        self.isE = 0  # is there an Earth?
        self.print_stuff = 0  # print information or not 
        self.read_spr()  # read in springs
        self.read_part() # read in  particles
        self.shape_vol = shape_vol   # shape model volume in m^3
        
    def read_spr(self):  # read springs from file
        i_index,j_index,ks,rs0,ggamma = read_springs(self.mdir,self.froot,self.index)
        self.i_index = i_index # spring connects particles i and j
        self.j_index = j_index
        self.rs0 = rs0  # is a vector, rest length 
        self.ks = ks  # is probably a single constant, spring constant
        self.ggamma = ggamma  # is probably a single constant, damping coefficient 
        
    def read_part(self): # read particles from file 
        xarr,yarr,zarr,vxarr,vyarr,vzarr,marr = read_particles(self.mdir,\
                                        self.froot,self.index)
        if (marr[-1] != marr[0]):
            self.isE=1
            self.ME = marr[-1] # EARTH
            self.posE = np.array([xarr[-1],yarr[-1],zarr[-1]])
            self.velE = np.array([vxarr[-1],vyarr[-1],vzarr[-1]])
            xarr = xarr[0:-1]
            yarr = yarr[0:-1]
            zarr = zarr[0:-1]
            marr = marr[0:-1]
            vxarr = vxarr[0:-1]
            vyarr = vyarr[0:-1]
            vzarr = vzarr[0:-1]
            if (self.print_stuff ==1):
                print('Earth is present')
            
        self.pos = np.zeros((len(xarr),3))
        self.vel = np.zeros((len(xarr),3))
        self.pos[:,0] = xarr
        self.pos[:,1] = yarr
        self.pos[:,2] = zarr
        self.vel[:,0] = vxarr
        self.vel[:,1] = vyarr
        self.vel[:,2] = vzarr
        self.m = marr[0] #  assuming all the same mass 
        
    # compute stress tensor at each node by hand
    def compute_stress_tensor(self):  
        n_p = self.pos.shape[0]  #number of particles 
        #print(n_p)
        ns = len(self.i_index)   # number of springs 
        sigma = np.zeros((n_p,3,3))
        for k in range(ns): # loop over springs
            pi = self.i_index[k]  # spring connects particle pi to pj
            pj = self.j_index[k]
            Force,Lvec = self.spring_force_one(k)
            for i in range(3):
                for j in range(3):
                    sigma[pi,i,j] += Force[i]*Lvec[j]  
                    sigma[pj,i,j] += Force[i]*Lvec[j]
        nodevol = self.shape_vol/n_p  #  not 100% sure about normalization!
        self.sigma = sigma/nodevol # now has correct units 
        if (self.print_stuff ==1):
            print('stress tensors computed')
        self.get_stress_eigs(); # compute the eigenvalues of the stress tensor
  
    # compute force and length vector for 1 spring k (used in computation of stress tensor)
    def spring_force_one(self,k):
        pi = self.i_index[k]
        pj = self.j_index[k]
        r_i = np.squeeze(self.pos[pi,:])
        r_j = np.squeeze(self.pos[pj,:])
        #print(r_i,r_j)
        dr = r_i - r_j  # is  a vector 
        L = np.sqrt(np.sum(dr*dr))  # current spring length 
        ks  = self.ks  # spring constant 
        rs0 = self.rs0[k]  # rest length 
        fac = -ks*(L-rs0)/L 
        Force = fac*dr 
        Lvec = dr 
    
        if (self.ggamma>0.0):   # damping force too!
            v_i = np.squeeze(self.vel[pi,:])
            v_j = np.squeeze(self.vel[pj,:])
            dv = v_i - v_j
            dLdt = np.sum(dr*dv)/L
                # divide dL/dt by L to get strain rate
            mbar = self.m/2; # reduced mass
            dampfac  = self.ggamma*mbar*dLdt/L;
               # factor L here to normalize dr
            Force -= dampfac*dr
    
        return Force, Lvec
    
    def rotate_to_sim(self,sim1): # rotate our simulation to the orientation of sim1
        q,rssd = R.align_vectors(sim1.pos, self.pos, weights=None)
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.align_vectors.html
        # scipy.spatial.transform.Rotation.align_vectors
        # compute rotation transformation q by comparing two sets of points
        # uses Kabsch transform 
        rotated_pos = q.apply(self.pos) # transform second sim to first one's orientation!
        self.pos = rotated_pos # update positions 
        rotated_vel = q.apply(self.vel)  # transform velocities to new orientation 
        self.vel = rotated_vel
        self.q = q # store rotation 
        if (self.isE ==1): # also rotate the Earth's position!
            rotated_posE= q.apply(self.posE)
            self.posE = rotated_posE
        if (self.print_stuff ==1):
            print('positions rotated')
        
    # get the eigenvalues of the stress tensor and their eigenvectors 
    def get_stress_eigs(self):
        n_p = self.pos.shape[0]  #number of particles 
        eigenvecs = self.sigma*0.0  # 3x3 matrix for each node 
        eigenvals = self.pos*0.0   # vector for each node
        for i in range(n_p):
            sig = np.squeeze(self.sigma[i,:,:]) # now 2d tensor
            w,v = linalg.eig(sig)  #compute eigenvalues and eigenvectors of stress tensor
            wsort,vsort = sort_eigs_3(w,v)  # sort the eigenvals and vecs big to little
            eigenvecs[i] = vsort
            eigenvals[i] = wsort 
        self.s_eigenvecs = eigenvecs
        self.s_eigenvals = eigenvals
        self.compute_J2()  # compute J2= invariant of deviatoric stress
        if (self.print_stuff ==1):
            print('eigenvectors and eigenvalues of stress tensor computed')
        
    # compute sqrt J2 = invariant of deviatoric stress
    def compute_J2(self):   # after computing eigenvalues of stress tensor 
        s1 = np.squeeze(self.s_eigenvals[:,0])  
        s2 = np.squeeze(self.s_eigenvals[:,1])
        s3 = np.squeeze(self.s_eigenvals[:,2])
        self.sJ2 = np.sqrt(((s1-s2)**2 + (s1-s3)**2 + (s2-s3)**2)/6)

    # compute moment of inertia tensor and its eigenvalues and eigenvectors 
    def compute_mom_inertia(self):
        I = np.zeros((3,3))
        n_p = self.pos.shape[0]  #number of particles
        for k in range(n_p):
            for i in range(3):   #diagonal part 
                r2 = self.pos[k,0]**2 + self.pos[k,1]**2 + self.pos[k,2]**2
                I[i,i] += r2
            for i in range(3):   # including off diagonal part 
                for j in range(3):
                    I[i,j] -= self.pos[k,i]*self.pos[k,j]
        I*= self.m
        self.I = I 
        w,v = linalg.eig(I)
        wsort,vsort = sort_eigs_3(w,v)  # sort the eigenvals and vecs big to little
        self.I_eigenvals = wsort
        self.I_eigenvecs = vsort
        #print('moment of inertia computed')
        
    def rotate_to_principal(self):  # rotate to principal axis coordinate system
        #  compute_mom_inertia 
        self.compute_mom_inertia()  # compute moment of inertia along with its eigenvectors
        #### xxxx there might be a mistake here! if I am not rotating to principal axes
        q =  R.from_matrix(self.I_eigenvecs) # create rotation from eigenvectors
        # does this rotation have det 1 necessarily,  maybe not.  Perhaps that's the problem
        # q is a rotation class thing in python
        rotated_pos = q.apply(self.pos) # do rotation!
        xpmax = np.max( rotated_pos[:,0])
        xmmax = np.max(-rotated_pos[:,0])
        if (xmmax > xpmax):  # we want to flip x direction 
            qflip = R.from_rotvec(np.pi * np.array([0, 0, 1.]))  # rotate about z by 180deg
            q = qflip*q  # apply flip to our rotation 
            rotated_pos = q.apply(self.pos) # reapply rotation, including flip!
            print('body flipped')
        self.pos = rotated_pos # update positions 
        rotated_vel = q.apply(self.vel)
        self.vel = rotated_vel
        self.q_I = q # store rotation 
        if (self.isE ==1): # also rotate the Earth's position!
            rotated_posE= q.apply(self.posE)
            self.posE = rotated_posE
        if (self.print_stuff ==1):
            print('rotated to principal axis frame')
        self.compute_mom_inertia() #  recompute moments of inertia tensor
        
    def compute_L(self):   # compute angular momentum and spin vector 
        self.L = np.zeros(3)
        n_p = self.pos.shape[0]  #number of particles
        for i in range(n_p):  # compute spin angular momentum 
            rvec = np.squeeze(self.pos[i])
            vvec = np.squeeze(self.vel[i])
            dl = np.cross(rvec,vvec)*self.m  # sum over m(r x v)
            self.L += dl; 
        self.compute_mom_inertia()  # compute moment of inertia 
        I_inv = linalg.inv(self.I)  # invert moment of inertia 
        self.Omega = np.matmul(I_inv,self.L)  # compute spin vector 
        
    # compute spring lengths
    def compute_spring_lengths(self):
        ns = len(self.i_index)
        self.Larr = np.zeros(ns)
        for k in range(ns):  #loop over springs 
            pi = self.i_index[k]
            pj = self.j_index[k]
            r_i = np.squeeze(self.pos[pi,:])
            r_j = np.squeeze(self.pos[pj,:])
            dr = r_i - r_j  # is  a vector 
            L = np.sqrt(np.sum(dr*dr))  # current spring length 
            self.Larr[k] = L
        


# globals for figures 
fig_lim = 213; # limit of x,y,z ranges to show in figures 
min_circle_ratio=0.04 # helps with concavity in triangulation 
Evec_mag = 150.0;  # size of vector pointing to Earth

# plot components of stress tensor
# arguments:
#  sim which simulation to use
#  sim1 is used to make a difference  if subit==1, else ignored
#  vmin,vmax ranges for xx yy zz components
#  vmin2,vmax2, ranges for xy, yz, xz components
#  zstep1,zstep2, steps for colorbars
#  ofile output file name for png image
#  tlabel  a time label in hours from peri 
def plt_sig6(sim,sim1,subit,vmin,vmax,zstep1,vmin2,vmax2,zstep2,tlabel,ofile):
    fig, axarr = plt.subplots(6,3,figsize=(4.3,7),sharex=True,sharey=True,dpi=100)
    plt.subplots_adjust(wspace=0,hspace=0)
    
    width = 20.   #could be adjusted choose a width of plane to contain points 
    #min_circle_ratio=0.03 # helps with tricontour concavity 
    
    if (sim.isE ==1):
        dE = np.sqrt(np.sum(sim.posE*sim.posE))
        Evec = sim.posE/dE * Evec_mag # to display Earth direction 
    
    pos = sim.pos
    if (subit==1):
        sigma= sim.sigma  - sim1.sigma  # subtract out first one 
        #cmap_e =  mpl.cm.bwr; cmap_e.set_under('blue'); cmap_e.set_over('red')
    else:
        sigma = sim.sigma #- sim1.sigma
        #cmap_e =  mpl.cm.pink; cmap_e.set_under('black'); cmap_e.set_over('white')
    
    for i in range(6):
        for j in range(3):
            axarr[i,j].set_aspect('equal')
    
    #fig_lim = 212;
    axarr[0,0].set_xlim([-fig_lim,fig_lim])
    axarr[0,0].set_ylim([-fig_lim,fig_lim])
    axarr[0,0].set_xticks([])
    axarr[0,0].set_yticks([])

    for k in range(6):

        if (k==0):
            i=0; j=0; sig_str = 'xx'
        if (k==1):
            i=1; j=1; sig_str = 'yy'
        if (k==2):
            i=2; j=2; sig_str = 'zz'
        if (k==3):
            i=0; j=1; sig_str = 'xy'
        if (k==4):
            i=0; j=2; sig_str = 'xz'
        if (k==5):
            i=1; j=2; sig_str = 'yz'
         
        cmap_e =  mpl.cm.bwr; cmap_e.set_under('blue'); cmap_e.set_over('red')
        if (k <3):
            levs = np.arange(vmin,vmax,zstep1)
            if (subit==0):
                cmap_e =  mpl.cm.pink; cmap_e.set_under('black'); cmap_e.set_over('white')
        else:
            levs = np.arange(vmin2,vmax2,zstep2)
            #cmap_e =  mpl.cm.bwr; cmap_e.set_under('blue'); cmap_e.set_over('red')
            
        xl=-156;yl=172
        axarr[k,0].text(xl,yl,sig_str,fontsize=10,ha='center',va='center')
            
        ypos = np.squeeze(pos[:,1])  
        ii = (np.abs(ypos)<width) # restrict to near y=0
        x = pos[ii,0]; y = pos[ii,1]; z = pos[ii,2]
        triang = tri.Triangulation(x, z)
        mask = tri.TriAnalyzer(triang).get_flat_tri_mask(min_circle_ratio)
        triang.set_mask(mask)
        axarr[k,0].tricontourf(triang,sigma[ii,i,j],levels=levs,cmap=cmap_e,extend='both')

        xpos = np.squeeze(pos[:,0])
        ii = (np.abs(xpos)<width) # restrict to near x=0
        x = pos[ii,0]; y = pos[ii,1]; z = pos[ii,2]
        triang = tri.Triangulation(y, z)
        mask = tri.TriAnalyzer(triang).get_flat_tri_mask(min_circle_ratio)
        triang.set_mask(mask)
        axarr[k,1].tricontourf(triang,sigma[ii,i,j],levels=levs,cmap=cmap_e,extend='both')
        
        zpos = np.squeeze(pos[:,2])
        ii = (np.abs(zpos)<width) # restrict to near z=0
        x = pos[ii,0]; y = pos[ii,1]; z = pos[ii,2]
        triang = tri.Triangulation(y, x)
        mask = tri.TriAnalyzer(triang).get_flat_tri_mask(min_circle_ratio)
        triang.set_mask(mask)
        im=axarr[k,2].tricontourf(triang,sigma[ii,i,j],levels=levs,cmap=cmap_e,extend='both')
        
        if (k==0):
            cbar = plt.colorbar(im,ax=axarr[0:3,:],shrink=0.9,aspect=30,label='Stress (Pa)')
        if (k==3):   
            cbar = plt.colorbar(im,ax=axarr[3:6,:],shrink=0.9,aspect=30,label='Stress (Pa)')
    
    if (sim.isE ==1):
        for k in range(6): # show a segment pointing toward Earth
            axarr[k,0].plot([0,Evec[0]],[0,Evec[2]],'c-',lw=1)  # point toward the Earth!
            axarr[k,1].plot([0,Evec[1]],[0,Evec[2]],'c-',lw=1)
            axarr[k,2].plot([0,Evec[1]],[0,Evec[0]],'c-',lw=1)
    
    axarr[5,0].set_ylabel('z')
    axarr[5,0].set_xlabel('x')
    axarr[5,1].set_xlabel('y')
    axarr[5,2].set_xlabel('y')
    axarr[5,2].text(fig_lim,0,'x',rotation='vertical')
    #ax_r = axarr[5,2].twinx() bombs
    #ax_r.set_ylabel('y')
    
    x0=-40; y0=-170.0; x1 = x0+200.0;
    axarr[0,0].plot([x0,x1],[y0,y0],'k-',lw=2)   # length bar 
    axarr[0,0].text((x0+x1)/2,y0,'200 m',ha='center',va='bottom',fontsize=8)
    # show scale bar 
    
    tlabel_full = r'$t-t_{peri}$=' + tlabel + ' hour'
    axarr[0,1].text(0,240,tlabel_full,ha='center',va='center')
    
    if (len(ofile)>3):
        plt.savefig(ofile,dpi=200)

# plot deviatoric stress 
def plt_J2(sim,sim1,subit,vmin,vmax,zstep,tlabel,ofile):
    fig, axarr = plt.subplots(1,3,figsize=(6,2),sharex=True,sharey=True)
    plt.subplots_adjust(wspace=0,hspace=0,left=0.04,right=0.96,top=0.99,bottom=0.04)
    pos= sim.pos 
    if (subit==1):
        sJ2 = sim.sJ2 - sim1.sJ2
        #cmap_b =  mpl.cm.BrBG; cmap_b.set_under('brown'); cmap_b.set_over('green')
        # create nice colormap
        top = mpl.colormaps['BuGn_r'].resampled(128)
        bottom = mpl.colormaps['YlOrRd'].resampled(128)
        newcolors = np.vstack((top(np.linspace(0, 1, 128)),\
                       bottom(np.linspace(0, 1, 128))))
        cmap_b = mpl.colors.ListedColormap(newcolors, name='BuGnYlOrRd')
        cmap_b.set_under('green'); cmap_b.set_over('red')  # extend it
    else:
        sJ2 = np.sqrt(sim.sJ2)
        cmap_b =  mpl.cm.hot;  cmap_b.set_under('black'); cmap_b.set_over('white')
        # extended 
        
    #fig_lim = 212;
    axarr[0].set_xlim([-fig_lim,fig_lim])
    axarr[0].set_ylim([-fig_lim,fig_lim])
    axarr[0].set_xticks([])
    axarr[0].set_yticks([])
    
    width = 20.   # width of plane containing points that are plotted
    axarr[0].set_aspect('equal')
    axarr[1].set_aspect('equal')
    axarr[2].set_aspect('equal')
    
    levs  = np.arange(vmin,vmax,zstep)
    
    if (sim.isE ==1):
        dE = np.sqrt(np.sum(sim.posE*sim.posE)) # Earth!
        Evec = sim.posE/dE * Evec_mag;
    
    ypos = np.squeeze(pos[:,1])
    ii = (np.abs(ypos)<width)
    x = np.squeeze(pos[ii,0]); y = np.squeeze(pos[ii,1]); z = np.squeeze(pos[ii,2])
    triang = tri.Triangulation(x, z)
    mask = tri.TriAnalyzer(triang).get_flat_tri_mask(min_circle_ratio)
    triang.set_mask(mask)
    axarr[0].tricontourf(triang,sJ2[ii],levels=levs,cmap=cmap_b,extend='both')
    
    xpos = np.squeeze(pos[:,0])
    ii = (np.abs(xpos)<width)
    x = np.squeeze(pos[ii,0]); y = np.squeeze(pos[ii,1]); z = np.squeeze(pos[ii,2])
    triang = tri.Triangulation(y, z)
    mask = tri.TriAnalyzer(triang).get_flat_tri_mask(min_circle_ratio)
    triang.set_mask(mask)
    axarr[1].tricontourf(triang,sJ2[ii],levels=levs,cmap=cmap_b,extend='both')
    
    zpos = np.squeeze(pos[:,2])
    ii = (np.abs(zpos)<width)
    x = np.squeeze(pos[ii,0]); y = np.squeeze(pos[ii,1]); z = np.squeeze(pos[ii,2])
    triang = tri.Triangulation(y, x)
    mask = tri.TriAnalyzer(triang).get_flat_tri_mask(min_circle_ratio)
    triang.set_mask(mask)
    im= axarr[2].tricontourf(triang,sJ2[ii],levels=levs,cmap=cmap_b,extend='both')

    if (sim.isE ==1):
        axarr[0].plot([0,Evec[0]],[0,Evec[2]],'c-',lw=1) # Earth
        axarr[1].plot([0,Evec[1]],[0,Evec[2]],'c-',lw=1)
        axarr[2].plot([0,Evec[1]],[0,Evec[0]],'c-',lw=1)

    cstring = r'$\sqrt{J_2}$(Pa)'
    if (subit==1):
        cstring = r'$\Delta$' + cstring
        
    plt.colorbar(im,ax=axarr,shrink=0.8,label=cstring)
    
    x0=-40; y0=-180.0; x1 = x0+200.0;
    axarr[0].plot([x0,x1],[y0,y0],'k-',lw=2)  # put a 200 m bar on the plot
    axarr[0].text((x0+x1)/2,y0,'200 m',ha='center',va='bottom',fontsize=8)
    axarr[0].set_ylabel('z')
    axarr[0].set_xlabel('x')
    axarr[1].set_xlabel('y')
    axarr[2].set_xlabel('y')
    axarr[2].text(fig_lim,0,'x',rotation='vertical')
    
    tlabel_full = r'$t-t_{peri}$=' + tlabel + ' hour'
    axarr[0].text(0,230,tlabel_full,ha='center',va='center')
    
    if (len(ofile)>3):
        plt.savefig(ofile,dpi=300)


# plot deformation with a quiver plot!
# arguments:
#    Qscale: you can adjust to get arrow the right size
#    Qkey_len:  for quiverkey
#    sim, sim1: simulations. Show positions computed from sim-sim1 
#    tlabel is a time lable
#    ofile, string for png if you want a figure output
def plt_def3(sim,sim1,Qscale,Qkey_len,tlabel,ofile):
    
    fig, axarr = plt.subplots(1,3,figsize=(6.0,2.5),sharex=True,sharey=True,dpi=200)
    plt.subplots_adjust(wspace=0,hspace=0,left=0.04,right=0.96,top=0.97,bottom=0.08)

    for i in range(3):
        axarr[i].set_aspect('equal')
    #fig_lim = 212;
    axarr[0].set_xlim([-fig_lim,fig_lim])
    axarr[0].set_ylim([-fig_lim,fig_lim])
    axarr[0].set_xticks([])
    axarr[0].set_yticks([])

    width = 20.   #could be adjusted choose a width of plane to contain points
    #min_circle_ratio=0.04 # helps with tricontour concavity

    pos = sim.pos; xpos = pos[:,0]; ypos = pos[:,1]; zpos = pos[:,2]
    dpos= sim.pos  - sim1.pos  # subtract out other simulation
    
    if (sim.isE ==1):
        dE = np.sqrt(np.sum(sim.posE*sim.posE)) # Earth!
        Evec = sim.posE/dE * Evec_mag;
    
    # Quiver props
    #Qscale = 5e-4; 
    Qwidth=2; Qkeyx = -150; Qkeyy=200; 
    #Qkey_len=0.03; 
    Qkey_lab = '{:.1f}cm'.format(Qkey_len*100)
    
    ii = (np.abs(ypos)<width) # restrict to near y=0
    x  = np.squeeze(pos[ii,0]) ; y  = np.squeeze(pos[ii,1]);   z = np.squeeze(pos[ii,2])
    dx = np.squeeze(dpos[ii,0]); dy = np.squeeze(dpos[ii,1]); dz = np.squeeze(dpos[ii,2])
    Q = axarr[0].quiver(x,z,dx,dz,units='x',scale=Qscale,pivot='mid',width=Qwidth)
    qk = axarr[0].quiverkey(Q, Qkeyx, Qkeyy, Qkey_len, '', labelpos='S',\
                   coordinates='data')
    axarr[0].text(Qkeyx,Qkeyy-20,Qkey_lab,fontsize=8,ha='center',va='center')
    
    ii = (np.abs(xpos)<width) # restrict to near x=0
    x  = np.squeeze(pos[ii,0]) ; y  = np.squeeze(pos[ii,1]);   z = np.squeeze(pos[ii,2])
    dx = np.squeeze(dpos[ii,0]); dy = np.squeeze(dpos[ii,1]); dz = np.squeeze(dpos[ii,2])
    Q= axarr[1].quiver(y,z,dy,dz,units='x',scale=Qscale,pivot='mid',width=Qwidth)
    
    ii = (np.abs(zpos)<width) # restrict to near z=0
    x  = np.squeeze(pos[ii,0]) ; y  = np.squeeze(pos[ii,1]);   z = np.squeeze(pos[ii,2])
    dx = np.squeeze(dpos[ii,0]); dy = np.squeeze(dpos[ii,1]); dz = np.squeeze(dpos[ii,2])
    axarr[2].quiver(y,x,dy,dx,units='x',scale=Qscale,pivot='mid',width=Qwidth)

    if (sim.isE ==1):
        axarr[0].plot([0,Evec[0]],[0,Evec[2]],'c-',lw=1) # Earth
        axarr[1].plot([0,Evec[1]],[0,Evec[2]],'c-',lw=1) # Earth
        axarr[2].plot([0,Evec[1]],[0,Evec[0]],'c-',lw=1) # Earth
    
    x0=-40; y0=-190.0; x1 = x0+200.0;
    axarr[0].plot([x0,x1],[y0,y0],'k-',lw=2)  # put a 200 m bar on the plot
    axarr[0].text((x0+x1)/2,y0,'200 m',ha='center',va='bottom',fontsize=8)
    axarr[0].set_ylabel('z')
    axarr[0].set_xlabel('x')
    axarr[1].set_xlabel('y')
    axarr[2].set_xlabel('y')
    axarr[2].text(fig_lim,0,'x',rotation='vertical')
    
    tlabel_full = r'$t-t_{peri}$=' + tlabel + ' hour'
    axarr[0].text(0,230,tlabel_full,ha='center',va='center')

    if (len(ofile)>3):
        plt.savefig(ofile,dpi=300)


