# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 17:55:54 2011

@author: Sat Kumar Tomer
@website: http://www.ambhas.com/tools/html/krige_8py_source.html
@email: satkumartomer@gmail.com

"""

# import required modules
import numpy as np
import matplotlib.pylab as plt


class OK:
    """
    This performs the ordinary kriging
    Input:
        x: x vector of location
        Y: y vector of location
        z: data vector at location (x,y)

    Output:
        None

    Methods:
        variogram: estimate the variogram

    """
    def __init__(self,x,y,z):
        self.x = x.flatten()
        self.y = y.flatten()
        self.z = z.flatten()

    def variogram(self, var_type='averaged', n_lag=9):
        """
        var_type: averaged or scattered
        """

        x = self.x
        y = self.y
        z = self.z
        # make the meshgrid
        X1,X2 = np.meshgrid(x,x)
        Y1,Y2 = np.meshgrid(y,y)
        Z1,Z2 = np.meshgrid(z,z)

        D = np.sqrt((X1 - X2)**2 + (Y1 - Y2)**2)

        G = 0.5*(Z1 - Z2)**2
        indx = range(len(z))
        C,R = np.meshgrid(indx,indx)
        G = G[R>C]

        self.D = D
        DI = D[R > C]

        # group the variogram
        # the group are formed based on the equal number of bin
        total_n = len(DI)
        group_n = int(total_n/n_lag)
        sor_i = np.argsort(DI)[::-1]

        DE = np.empty(n_lag)
        GE = np.empty(n_lag)
        for i in range(n_lag):
            if i<n_lag-1:
                DE[i] = DI[sor_i[group_n*i:group_n*(i+1)]].mean()
                GE[i] = G[sor_i[group_n*i:group_n*(i+1)]].mean()

            else:
                DE[i] = DI[sor_i[group_n*i:]].mean()
                GE[i] = G[sor_i[group_n*i:]].mean()

        if var_type == 'scattered':
            return DI,G
        elif var_type == 'averaged':
            return DE,GE
        else:
            raise ValueError('var_type should be either averaged or r')

    def vario_model(self, lags, model_par, model_type='linear'):
        """
        Input:
            model_type : the type of variogram model
                             spherical
                             linear
                             exponential
            model_par:  parameters of variogram model
                        this should be a dictionary
                        e.g. for shperical and exponential
                            model_par = {'nugget':0, 'range':1, 'sill':1}
                        for linear
                            model_par = {'nugget':0, 'slope':1}
        Output:
            G:  The fitted variogram model
        """

        if model_type == 'spherical':
            n = model_par['nugget']
            r = model_par['range']
            s = model_par['sill']
            l = lags
            G = n + (s*(1.5*l/r - 0.5*(l/r)**3)*(l<=r) + s*(l>r))

        elif model_type == 'linear':
            n = model_par['nugget']
            s = model_par['slope']
            l = lags
            G = n + s*l

        elif model_type == 'exponential':
            n = model_par['nugget']
            r = model_par['range']
            s = model_par['sill']
            l = lags
            G = n + s*(1 - np.exp(-3*l/r))

        else:
            raise ValueError('model_type should be spherical or linear or ntial')

        return G

    def int_vario(self, Xg, Yg, model_par, model_type):
        """
        this computes the integral of the variogram over a square
        using the Monte Carlo integration method

        this works only for two dimensional grid

        Input:
            Xg:     x location where krigged data is required
            Yg:     y location whre kirgged data is required
            model_par: see the vario_model
            model_type: see the vario_model
        """
        avg_vario = np.empty((len(self.x), (len(Xg)-1)*(len(Yg)-1)))
        for k in range(len(self.x)):

            avg_vario_ens = np.empty((len(Xg)-1, len(Yg)-1))
            for i in range(len(Xg)-1):
                for j in range(len(Yg)-1):
                    Xg_rand = Xg[i]+np.random.rand(10)*(Xg[i+1]-Xg[i])
                    Yg_rand = Yg[j]+np.random.rand(10)*(Yg[j+1]-Yg[j])

                    DOR = ((self.x[k] - Xg_rand)**2 + (self.y[k] -
d)**2)**0.5
                    avg_vario_ens[i,j] = self.vario_model(DOR, model_par,
type).mean()
            avg_vario[k,:] = avg_vario_ens.flatten()
        return avg_vario

    def krige(self, Xg, Yg, model_par, model_type):
        """
        Input:
            Xg:     x location where krigged data is required
            Yg:     y location whre kirgged data is required
            model_par: see the vario_model
            model_type: see the vario_model

        Attributes:
            self.Zg : krigged data
            self.s2_k = variance in the data

        """

        # set up the Gmod matrix
        n = len(self.x)
        Gmod = np.empty((n+1,n+1))
        Gmod[:n, :n] = self.vario_model(self.D, model_par, model_type)

        Gmod[:,n] = 1
        Gmod[n,:] = 1
        Gmod[n,n] = 0

        Gmod = np.matrix(Gmod)

        # inverse of Gmod
        Ginv = Gmod.I

        Xg = Xg.flatten()
        Yg = Yg.flatten()
        Zg = np.empty(Xg.shape)
        s2_k = np.empty(Xg.shape)

        for k in range(len(Xg)):

            DOR = ((self.x - Xg[k])**2 + (self.y - Yg[k])**2)**0.5
            GR = np.empty((n+1,1))

            GR[:n,0] = self.vario_model(DOR, model_par, model_type)

            GR[n,0] = 1
            E = np.array(Ginv * GR )
            Zg[k] = np.sum(E[:n,0]*self.z)
            s2_k[k] = np.sum(E[:n,0]*GR[:n,0])+ E[n, 0]

        self.Zg = Zg
        self.s2_k = s2_k

    def block_krige(self, Xg, Yg, model_par, model_type):
        """
        Input:
            Xg:     x location where krigged data is required
            Yg:     y location whre krigged data is required
            model_par: see the vario_model
            model_type: see the vario_model

        Attributes:
            self.Zg : krigged data
            self.s2_k = variance in the data

        """

        # set up the Gmod matrix
        n = len(self.x)
        Gmod = np.empty((n+1,n+1))
        Gmod[:n, :n] = self.vario_model(self.D, model_par, model_type)

        Gmod[:,n] = 1
        Gmod[n,:] = 1
        Gmod[n,n] = 0

        Gmod = np.matrix(Gmod)

        # inverse of Gmod
        Ginv = Gmod.I

        Xg = Xg.flatten()
        Yg = Yg.flatten()


        avg_vario = self.int_vario(Xg, Yg, model_par, model_type)
        Zg = np.empty(avg_vario.shape[1])
        s2_k = np.empty(avg_vario.shape[1])

        for k in range(avg_vario.shape[1]):

            GR = np.empty((n+1,1))
            GR[:n,0] = avg_vario[:,k]
            GR[n,0] = 1
            E = np.array(Ginv * GR )
            Zg[k] = np.sum(E[:n,0]*self.z)
            s2_k[k] = np.sum(E[:n,0]*GR[:n,0])+ E[n, 0]

        self.Zg = Zg.reshape(len(Xg)-1, len(Yg)-1)
        self.s2_k = s2_k.reshape(len(Xg)-1, len(Yg)-1)

if __name__ == "__main__":
    # generate some sythetic data
    x = np.random.rand(20)
    y = np.random.rand(20)
    z = 0.0*np.random.normal(size=20)+x+y

    foo = OK(x,y,z)
    #ax,ay = foo.variogram('scattered')
    ax,ay = foo.variogram()

    plt.plot(ax,ay,'ro')

    lags = np.linspace(0,5)
    model_par = {}
    model_par['nugget'] = 0
    model_par['range'] = 1
    model_par['sill'] = 2.0

    G = foo.vario_model(lags, model_par, model_type = 'exponential')
    plt.plot(lags, G, 'k')
    plt.show()

    Rx = np.linspace(-1,1,1050)
    Ry = np.linspace(0,1,750)
    XI,YI = np.meshgrid(Rx,Ry)
    foo.krige(XI, YI, model_par, 'exponential')

    plt.matshow(foo.Zg.reshape(750,1050))
    plt.show()

#    # block kriging
#    xg = np.linspace(0,1,5)
#    yg = np.linspace(0,1,8)
#    foo.block_krige(xg, yg, model_par, model_type = 'exponential')
#    plt.imshow(foo.s2_k, extent=(0,1,0,1))
#    plt.imshow(foo.Zg, extent=(0,1,0,1))
#    plt.matshow(foo.Zg)
#    plt.matshow(foo.s2_k)
#    plt.colorbar()
#    plt.plot(x,y, 'ro')
#    plt.show()
