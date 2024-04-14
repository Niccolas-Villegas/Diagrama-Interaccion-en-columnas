import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from units import *



class rebar:
    def __init__(self, B):
        self.B = B
    
    def diametro(self, db):
        '''
        Determina el diametro del acero longitudinal segun su numero de barra
        '''
        if db == 2:
            db = 6.4*mm
        elif db == 3:
            db = 9.5*mm
        elif db == 4:
            db = 12.7*mm
        elif db == 5:
            db = 15.9*mm
        elif db == 6:
            db = 19.1*mm
        elif db == 8:
            db = 25.4*mm
        else:
            db = 0.00
        return db

    def area(self, db):
        '''
        Determina el area del acero longitudinal segun su numero de barra
        '''
        if db == 2:
            As = 0.32*cm**2
        elif db == 3:
            As = 0.71*cm**2
        elif db == 4:
            As = 1.27*cm**2
        elif db == 5:
            As = 1.98*cm**2
        elif db == 6:
            As = 2.84*cm**2
        elif db == 8:
            As = 5.10*cm**2
        else:
            As = 0.00
        return As
    
    def As(self):
        '''
        Determina el area del acero de refuerzo por cada tramo de barras segun la direccion analizada
        '''
        nbx, nby = self.B.shape
        A = np.zeros((nbx, nby))
        for i in range(nbx):
            for j in range(nby):
                A[i][j] = self.area(self.B[i][j])
        return A

class variables_generales:
    # Segun ACI 318-19 - Tabla 22.2.2.4.3.
    def fs_valor(self, es):
        Es = (2*10**6)*kgf/cm**2
        if es > 0.0021:
            fs = 4200*kgf/cm**2
        elif es < -0.0021:
            fs = -4200*kgf/cm**2
        else:
            fs = Es*es
        return fs
    
    # Segun ACI 318-19 - Tabla 22.2.2.4.3.
    def beta(self, fc):
        if 175*kgf/cm**2 <= fc <= 280*kgf/cm**2:
            beta = 0.85
        elif 280*kgf/cm**2 < fc < 550*kgf/cm**2:
            beta = 0.85 - 0.05*(self.fc-280*kgf/cm**2)/70*kgf/cm**2
        elif fc > 550*kgf/cm**2:
            beta = 0.65
        return beta

class Diagrama_Interaccion_Rect:
    def __init__(self, b, h, rec, d_est, fc, fy, Rebar, P_dem, Mx_dem, My_dem):
        self.b = b
        self.h = h
        self.bo = b
        self.ho = h
        self.rec = rec
        self.d_est = d_est
        self.fc = fc
        self.fy = fy
        self.Rebar = Rebar
        self.P_dem = P_dem
        self.Mx_dem = Mx_dem
        self.My_dem = My_dem
        
        for eje in ['X','Y']:
            self.ubicacion, self.As = self.direccion(eje)
            if eje == 'X':
                self.y_loc = self.ubicacion
                self.Pnx, self.Mnx, self.ACI_Pnx, self.ACI_Mnx, self.E060_Pnx, self.E060_Mnx= self.curva_reducida()
                self.Pbx, self.Mbx = self.zona_balance()
            else:
                self.x_loc = self.ubicacion
                self.Pny, self.Mny ,self.ACI_Pny, self.ACI_Mny , self.E060_Pny, self.E060_Mny= self.curva_reducida()
                self.Pby, self.Mby = self.zona_balance()
    
    def direccion(self, eje):
        '''
        Reorienta los inputs para obtener el diagrama de interaccion en la direccion a analizar
        '''            
        if eje == 'X':
            barras = np.max(self.Rebar, axis=1)
            As = np.sum(rebar(self.Rebar).As(), axis=1)
            n_capas = self.Rebar.shape[0]
        elif eje == 'Y':
            barras = np.max(self.Rebar, axis=0)
            As = np.sum(rebar(self.Rebar).As(), axis=0)
            self.b, self.h = self.h, self.b
            n_capas = self.Rebar.shape[1]
        d_barras = [rebar(self.Rebar).diametro(db) for db in barras]
        valor_fijo = self.h - 2*self.rec - 2*self.d_est
        espacio = (valor_fijo - np.sum(d_barras))/(n_capas-1)
        ubicacion = []
        for i in range(len(d_barras)):
            sum_barras = np.sum(d_barras[:i]) + d_barras[i]/2
            y = self.rec + self.d_est + sum_barras + i*espacio
            ubicacion.append(y)
        return ubicacion, As
    
    # Factor de reduccion segun ACI 318-19 - Tabla 21.2.2.
    def phi_ACI(self, es):
        '''
        Determina el factor de reduccion segun la norma ACI 318-19
        es: deformacion unitaria del refuerzo mas cercano a la cara en traccion 
        '''
        ety = 0.0021
        if abs(es) <= 0.0021:
            phi = 0.65
        elif 0.0021 < abs(es) <= 0.005:
            phi = 0.65 + 0.25*(abs(es)-ety)/(0.005-ety)
        else:
            phi = 0.90
        return phi
    
    # Factor de reduccion segun E.060 - Articulo 9.3.2.2.
    def phi_E060(self, Pn):
        '''
        Determina el factor de reduccion segun la norma E060
        Nota: Se aproximo el peralte efectivo como la altura menos 6cm
        Pn: Carga axial que se encuentra en Pmin y 0
        '''
        ety = 0.0021
        cb = (self.h-6*cm)/(1 + ety/0.003)
        Pnb = self.curva_nominal(cb)[0]
        Pmin = min((0.1*self.fc*self.b*self.h)/(0.7*tnf), Pnb)
        if Pn > Pmin:
            phi = 0.70
        elif Pn < 0:
            phi = 0.90
        else:
            phi = 0.90 - (0.90-0.70)*Pn/Pmin
        return phi
    
    def limites(self):
        '''
        Define los limites a compresion y traccion del diagrama de interaccion
        '''
        Pon = (0.85*self.fc*(self.b*self.h-sum(self.As)) + self.fy*sum(self.As))
        Tn = -self.fy*sum(self.As)
        return Pon/tnf, Tn/tnf

    def curva_nominal(self, c):
        beta = variables_generales().beta(self.fc)
        a = beta*c
        Fc = 0.85*self.fc*a*self.b
        sumFs, sumMs = 2*[0]
        for i in range(len(self.ubicacion)):
            es = 0.003*(c - self.ubicacion[i])/c
            fs = variables_generales().fs_valor(es)
            Fs = self.As[i]*fs
            sumFs = sumFs + Fs
            sumMs = sumMs + Fs*(self.h/2 - self.ubicacion[i])
        Pn = Fc + sumFs
        Mn = Fc*(self.h - a)/2 + sumMs
        return Pn/tnf, Mn/(tnf*m)

    def curva_ACI(self, c, Pn, Mn):
        '''
        Multiplica Pn y Mn por el factor phi según normativa ACI 318-19
        '''
        et = 0.003*(c - max(self.ubicacion))/c
        phi = self.phi_ACI(et)
        return phi*Pn, phi*Mn
    
    def curva_E060(self, Pn, Mn):
        '''
        Multiplica Pn y Mn por el factor phi según normativa E060
        '''
        phi = self.phi_E060(Pn)
        return phi*Pn, phi*Mn

    def curva_reducida(self):
        # Se inicia la lista de curvas con estos valores para cerrar los diagramas de interaccion
        Pon, Tn = self.limites()
        Pn_max = 0.80*Pon
        ACI_Pn_max = 0.65*Pn_max
        E060_Pn_max = 0.70*Pn_max
        Pnt = [Tn] + 2*[Pn_max]
        Mnt = 2*[0]
        ACI_Pnt = [0.90*Tn] + 2*[ACI_Pn_max]
        ACI_Mnt = 2*[0] 
        E060_Pnt = [0.90*Tn] + 2*[E060_Pn_max]
        E060_Mnt = 2*[0]
        
        # Obtencion de los diagramas de interaccion nominal y reducido (segun normativas)
        i = 1
        for c in np.linspace(1*cm, self.h, 30):
            Pn, Mn = self.curva_nominal(c)
            if Pn < Pn_max:
                Pnt.insert(i, Pn)
                Mnt.insert(i ,Mn)
                ACI_Pn, ACI_Mn = self.curva_ACI(c, Pn, Mn)
                ACI_Pnt.insert(i, ACI_Pn)
                ACI_Mnt.insert(i, ACI_Mn)
                E060_Pn, E060_Mn = self.curva_E060(Pn, Mn)
                E060_Pnt.insert(i, E060_Pn)
                E060_Mnt.insert(i, E060_Mn)
                i += 1
        
        # Interpolacion lineal para cerrar los diagramas de interaccion en Pn_max
        Mnf = Mnt[i-1] + (Mnt[i-1]-Mnt[i-2])/(Pnt[i-1]-Pnt[i-2])*(Pn_max-Pnt[i-1])
        Mnt.insert(i ,Mnf)
        
        ACI_Mnt.insert(i, 0.65*Mnf)
        E060_Mnt.insert(i, 0.70*Mnf)
        
        return Pnt, Mnt, ACI_Pnt, ACI_Mnt, E060_Pnt, E060_Mnt
    
    def zona_balance(self):
        '''
        Determina la zona de trancision o falla balanceada de la seccion
        Nota: Se aproximo el peralte efectivo como la altura menos 6cm
        '''
        ety = 0.0021
        cb = (self.h-6*cm)/(1 + ety/0.003)
        Pnb, Mnb = self.curva_nominal(cb)
        Pb = [Pnb, 0, Pnb]
        Mb = [-Mnb, 0, Mnb]
        
        return Pb, Mb
    
    def Ploteos(self, showACI : bool, showE060 : bool):
        Pnx = self.Pnx + self.Pnx[::-1]
        Pny = self.Pny + self.Pny[::-1]
        Mnx = self.Mnx + [-x for x in self.Mnx][::-1]
        Mny = self.Mny + [-y for y in self.Mny][::-1]
        ACI_Pnx = self.ACI_Pnx + self.ACI_Pnx[::-1]
        ACI_Pny = self.ACI_Pny + self.ACI_Pny[::-1]
        ACI_Mnx = self.ACI_Mnx + [-x for x in self.ACI_Mnx][::-1]
        ACI_Mny = self.ACI_Mny + [-y for y in self.ACI_Mny][::-1]
        E060_Pnx = self.E060_Pnx + self.E060_Pnx[::-1]
        E060_Pny = self.E060_Pny + self.E060_Pny[::-1]
        E060_Mnx = self.E060_Mnx + [-x for x in self.E060_Mnx][::-1]
        E060_Mny = self.E060_Mny + [-y for y in self.E060_Mny][::-1]
        
        # # Extracción de resultados
        # data1 = {'Pnx':Pnx, 'Mnx':Mnx, 'phi_Pnx':ACI_Pnx, 'phi_Mnx':ACI_Mnx}
        # df1 = pd.DataFrame(data1)
        # df1.to_csv('Resultados_x.csv')
        # data2 = {'Pny':Pny, 'Mny':Mny, 'phi_Pny':ACI_Pny, 'phi_Mny':ACI_Mny}
        # df2 = pd.DataFrame(data2)
        # df2.to_csv('Resultados_y.csv')

        fig, axs = plt.subplots(1, 2, figsize =(20,12))
        plt.rcParams['font.family'] = 'Segoe UI'    

        # Grafico de diagrama de interaccion en eje X
        legend_x = []
        axs[0].plot(Mnx, Pnx, 'k')
        legend_x.append('Capacidad Nominal')
        if showACI == True:
            axs[0].plot(ACI_Mnx, ACI_Pnx, 'b')
            legend_x.append('Capacidad ACI 318-19')
        if showE060 == True:
            axs[0].plot(E060_Mnx, E060_Pnx, 'r')
            legend_x.append('Capacidad E060')
        axs[0].plot(self.Mbx, self.Pbx, linestyle='--', color='silver') 
        for i in range(len(self.P_dem)):
            axs[0].plot(self.Mx_dem[i], self.P_dem[i], 'gx', markersize = 12)
        legend_x.extend(['Falla balanceada', 'Demanda'])
        axs[0].legend(legend_x, loc='lower right')
        axs[0].set_title('Diagrama de Interaccion - Eje X', size = 15)
        axs[0].grid(color = 'silver')
        axs[0].set_xlabel('Momento $ton.m$', size = 10)
        axs[0].set_ylabel('Carga Axial $ton$', size = 10)

        # Grafico de diagrama de interaccion en eje Y
        legend_y = []
        axs[1].plot(Mny, Pny, 'k')
        legend_y.append('Capacidad Nominal')
        if showACI == True:
            axs[1].plot(ACI_Mny, ACI_Pny, 'b')
            legend_y.append('Capacidad ACI 318-19')
        if showE060 == True:
            axs[1].plot(E060_Mny, E060_Pny, 'r')
            legend_y.append('Capacidad E060')
        axs[1].plot(self.Mby, self.Pby, linestyle='--', color='silver') 
        for i in range(len(self.P_dem)):
            axs[1].plot(self.My_dem[i], self.P_dem[i], 'gx', markersize = 12)
        legend_y.extend(['Falla balanceada', 'Demanda'])
        axs[1].legend(legend_y, loc='lower right')
        axs[1].set_title('Diagrama de Interaccion - Eje Y', size = 15)
        axs[1].grid(color = 'silver')
        axs[1].set_xlabel('Momento $ton.m$', size = 10)
        axs[1].set_ylabel('Carga Axial $ton$', size = 10)

        return plt.show()