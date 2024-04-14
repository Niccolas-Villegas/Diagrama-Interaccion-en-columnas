from P1_Diagrama_Interaccion import *
from units import *

# EJEMPLO DE APLICACION

# Geometría de la sección
# -----------------------
b = 30*cm               # base
h = 60*cm               # altura
rec = 4*cm              # recubrimiento
d_est = 9.5*mm          # diametro de estribo

# Propiedades del material
# ------------------------
fc = 280*kgf/cm**2      # Resistencia compresión concreto
fy = 4200*kgf/cm**2     # Resistencia fluencia acero

# Distribución del acero de refuerzo
# ----------------------------------
Rebar = np.array([[8, 8, 8, 8],
                  [8, 0, 0, 8],
                  [8, 8, 8, 8]])

# Solicitaciones de carga (tnf)
# -----------------------------
P_dem = [150, -100, 50, 0]       # Carga axial
Mx_dem = [20, -15.5, -36, 30]    # Momento flector X-X
My_dem = [10, -2, 8, -15]        # Momento flector Y-Y


x = Diagrama_Interaccion_Rect(b, h, rec, d_est, fc, fy, Rebar, P_dem, Mx_dem, My_dem)

x.Ploteos(showACI = True, showE060 = True)