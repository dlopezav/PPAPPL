

###################################### Frame classes to add in main Notebook  #########################################

import tkinter as tk
from tkinter import ttk, filedialog
from ttkbootstrap.constants import *
from ttkbootstrap.scrolled import ScrolledFrame
import pandas as pd
import math



class Notebook_Model_param(ttk.Frame):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook, style='cyborg')

        self.frames_systems = []

        self.frames = []

        self.Labelframes = [ttk.LabelFrame(self, bootstyle='info',text="MODEL")]

        self.labels = []

        self.entrys = []

        self.texts = []

        self.float_vars = []

        self.texts_values= ["time_t:", 0,
                            "total_time:", 0.005 ,
                            "deltat2:",0.005,
                            "initial value of turbocharger speed [unit?]:", 7522,
                            "initial pressure ratio compressor:",2.4,
                            "initial pressure ratio turbine:", 3.8,
                            "Initial_speed: ", 12.15,
                            "model time of the simulation:", 0.035,]
        
        
        for i in range(len(self.texts_values)):
            if i%2!=0:
                self.float_vars.append(tk.DoubleVar(value=self.texts_values[i]))
            else:
                self.texts.append(self.texts_values[i])


        for count, text in enumerate(self.texts):
                self.frames.append(ttk.Frame(self.Labelframes[0]))
                self.labels.append(ttk.Label(self.frames[count], text=text))
                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))

                
        
        self.pack(fill=BOTH, expand=YES)

        for i in range(len(self.Labelframes)):
            self.Labelframes[i].pack(fill=BOTH, pady=10, padx=20)

        self.update_idletasks()

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)

    def return_values(self):
        N_cycle = int(self.float_vars(1) // self.float_vars(2) +1)

        return  N_cycle,[var.get() for var in self.float_vars]
    


################ First system propeller ################################

class Notebook_propeller(ttk.Frame):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook, style='cyborg')
        
        self.scrollframe = ScrolledFrame(self)
    
        self.frames_systems = []

        self.frames = []

        self.Labelframes = [ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Constants"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="To be determined")]

        self.labels = []

        self.entrys = []

        self.texts = []

        self.float_vars = []

        self.texts_values =["Ship displacement volume [m3]:", 112404,
                            "Density of sea water [kg/m^3]:", 1026,
                            "Mass of ship [kg]:", 260000000,
                            "Diameter [m]:", 10,
                            "Thrust:", 0.2,
                            "dVs_dt = c_1 * V_s^2 [kg/m]:", 26250,
                            "Number of propeller blades, zp:", 6,
                            "Disk area coefficient, AE/Ao:", 0.82,
                            "Pitch to diameter ratio, p/Dp:", 0.93846154,
                            "Ship wake fraction, w, which is considered constant \n taking values in the range from 0.20:", 0.3]
        
        for i in range(len(self.texts_values)):
            if i%2!=0:
                self.float_vars.append(tk.DoubleVar(value=self.texts_values[i]))
            else:
                self.texts.append(self.texts_values[i])

        for count, text in enumerate(self.texts):
            if count < 2:
                # For the ones outside in the first labelframe
                self.frames.append(ttk.Frame(self.Labelframes[0]))
                self.labels.append(ttk.Label(self.frames[count], text=text))
                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
            else:
                self.frames.append(ttk.Frame(self.Labelframes[1]))
                self.labels.append(ttk.Label(self.frames[count], text=text))
                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))


        # Wageningen
        self.wageningen_frame = ttk.Frame(self.scrollframe)
        self.wageningen_label = ttk.Label(self.wageningen_frame, text="Wageningen coefficients imports:")
        self.wageningen_file_entry = ttk.Entry(self.wageningen_frame, state='readonly')
        self.browse_button = ttk.Button(self.wageningen_frame, text="Browse", command=self.browse_wageningen_file, style='info')
        
        self.pack(fill=BOTH, expand=YES)
        self.scrollframe.pack(fill=BOTH, expand=YES)   
        self.Labelframes[0].pack(fill=BOTH, pady=10, padx=20)
        self.Labelframes[1].pack(fill=BOTH, pady=10, padx=20)
        self.update_idletasks()

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)


        # Position wageningen frame and label
        self.wageningen_frame.pack(side=TOP, pady=5)
        self.wageningen_label.pack(side=TOP, padx=5)
        self.browse_button.pack(side=LEFT, padx=5)
        self.wageningen_file_entry.pack(side=LEFT, padx=5)


    def browse_wageningen_file(self):

        self.file_path = filedialog.askopenfilename(filetypes=[("Excel files", "*.xlsx;*.xls")])
        self.wageningen_file_entry.config(state='normal')
        self.wageningen_file_entry.delete(0, tk.END)
        self.wageningen_file_entry.insert(0, self.file_path)
        self.wageningen_file_entry.config(state='readonly')

        

  
    def return_values(self):

        df = pd.read_excel(self.file_path)
        n1_nump = df['n1'].to_numpy()
        ct_nump = df['ct'].to_numpy()
        s1_nump = df['s'].to_numpy()
        t1_nump = df['t'].to_numpy()
        u1_nump = df['u'].to_numpy()
        v1_nump = df['v'].to_numpy()

        n2_nump = df['n2'].to_numpy()
        cq_nump = df['cq'].to_numpy()
        s2_nump = df['s2'].to_numpy()
        t2_nump = df['t2'].to_numpy()
        u2_nump = df['u2'].to_numpy()
        v2_nump = df['v2'].to_numpy()

        m_hydro = self.float_vars(0)*self.float_vars(1)

        return m_hydro, n1_nump, ct_nump, s1_nump, t1_nump, u1_nump, v1_nump, n2_nump, cq_nump, s2_nump, t2_nump, u2_nump, v2_nump, \
               [var.get() for var in self.float_vars]
    

################ END FIRST SYSTEM ######################################################

################## Second system PID controller #########################################

class Notebook_PIDcontroller(ttk.Frame):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook, style='cyborg')
    
        self.frames_systems = []

        self.frames = []

        self.Labelframes = [ttk.LabelFrame(self, bootstyle='info',text="Constants"),
                            ttk.LabelFrame(self, bootstyle='info',text="To be determined")]

        self.labels = []

        self.entrys = []

        self.texts = []

        self.float_vars = []

        self.texts_values =["Number of engine cylinders:", 12,
                            "Constant for PID controller	kp:", 0.003,
                            "Constant for PID controller	kd:", 0,
                            "Constant for PID controller	ki:", 0.0015,
                            "Original injeciton time (at the beginning of the simulation):",0.01427586]
        
        for i in range(len(self.texts_values)):
            if i%2!=0:
                self.float_vars.append(tk.DoubleVar(value=self.texts_values[i]))
            else:
                self.texts.append(self.texts_values[i])

        for count, text in enumerate(self.texts):
            if count < 1:
                # For the ones outside in the first labelframe
                self.frames.append(ttk.Frame(self.Labelframes[0]))
                self.labels.append(ttk.Label(self.frames[count], text=text))
                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
            else:
                self.frames.append(ttk.Frame(self.Labelframes[1]))
                self.labels.append(ttk.Label(self.frames[count], text=text))
                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))


        # Steadystate
        self.steadystate_frame = ttk.Frame(self)
        self.steadystate_label = ttk.Label(self.steadystate_frame, text="Steadystate imports:")
        self.steadystate_file_entry = ttk.Entry(self.steadystate_frame, state='readonly')
        self.browse_button = ttk.Button(self.steadystate_frame, text="Browse", command=self.browse_steadystate_file, style='info')
        
        self.pack(fill=BOTH, expand=YES)   
        self.Labelframes[0].pack(fill=BOTH, pady=10, padx=20)
        self.Labelframes[1].pack(fill=BOTH, pady=10, padx=20)
        self.update_idletasks()

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)
                

        # Position steadystate frame and label
        self.steadystate_frame.pack(side=TOP, pady=5)
        self.steadystate_label.pack(side=TOP, padx=5)
        self.browse_button.pack(side=LEFT, padx=5)
        self.steadystate_file_entry.pack(side=LEFT, padx=5)


    def browse_steadystate_file(self):

        self.file_path = filedialog.askopenfilename(filetypes=[("Excel files", "*.xlsx;*.xls")])
        self.steadystate_file_entry.config(state='normal')
        self.steadystate_file_entry.delete(0, tk.END)
        self.steadystate_file_entry.insert(0, self.file_path)
        self.steadystate_file_entry.config(state='readonly')

        

  
    def return_values(self):

        #order of speed 
        steadystate= pd.read_excel(self.file_path)
        Nord = 2*math.pi/60*steadystate['me_ne'].to_numpy()

        return  Nord, [var.get() for var in self.float_vars]


########### END SECOND SYSTEM #############################################################



############## Third system 0d cycle begin ################################################

class Notebook_0dCycle(ttk.Frame):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook, style='cyborg')
    
        self.scrollframe = ScrolledFrame(self)

        self.frames_systems = []

        self.frames = []

        self.Labelframes = [ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Constants"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="To be confirmed"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Wiebe parameters"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Parameters combustion p167"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Engine reference point")]

        self.labels = []

        self.entrys = []

        self.texts = []

        self.float_vars = []

        self.texts_values= ["Delta t:", 0.0001,
                            "Number of revolution per cycle:", 1,
                            "Clearance volume [cm3]:",164671,
                            "Lower heating value :", 42707,
                            "wall surface [units ?]:",1.102343986,
                            "Stroke:", 3.468,
                            "l: ", 3468,
                            "ar:", 1734,
                            "Perfect gas constant:",8.314462,

                            "Tw [K]:",300+273.15,
                            "Pressure at inlet valve closing [bar]:", 4.1,
                            "Temperature at the inlet valve closing [K]:", 307,
                            "Initial masse in the cylinder [g]:",7462,
                            "Compression ratio:", 15,
                            "reference temperature [K]:", 298.15,
                            "Stochiometric coefficient:", 14.8,
                            "R / molar mass of exhaust", 0.1341,
                            "R / mass of air", 0.287,
                            "R / mass of fuel gas", 0.0489086,
                            
                            "acomb:",0.42,
                            "bcomb:",0.49,
                            "am:",0.59,
                            "bm:",0.21,
                            "cm:",-0.96,
                            "delta_m:",0.27,

                            "a_id:",0.39,
                            "b_id:",0.105,
                            "c_id:",3.12,
                            "mWiebe:",0.1,

                            "Air fuel ratio for reference point:",2.435696566131349,
                            "ignition delay for reference point",1,
                            "Mass at inlet valve closing for reference point [g]:",7462,
                            "Engine speed for reference point [rpm]:",80,
                            "Wiebe curve form coefficient for reference point:",0.1,
                            "Combustion duration for reference point:",0.3667076110839844*60/(2*math.pi),
                            "SMFR:",14500,
                            "topen:",0.0005,
                            "tclose",0.0005]
        
        
        for i in range(len(self.texts_values)):
            if i%2!=0:
                self.float_vars.append(tk.DoubleVar(value=self.texts_values[i]))
            else:
                self.texts.append(self.texts_values[i])


        for count, text in enumerate(self.texts):
            if count < 9:
                # For the ones outside in the first labelframe
                self.frames.append(ttk.Frame(self.Labelframes[0]))
                self.labels.append(ttk.Label(self.frames[count], text=text))
                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
            else:
                if count<19:
                    self.frames.append(ttk.Frame(self.Labelframes[1]))
                    self.labels.append(ttk.Label(self.frames[count], text=text))
                    self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                else:
                    if count<25:
                        self.frames.append(ttk.Frame(self.Labelframes[2]))
                        self.labels.append(ttk.Label(self.frames[count], text=text))
                        self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                    else:
                        if count<29:
                            self.frames.append(ttk.Frame(self.Labelframes[3]))
                            self.labels.append(ttk.Label(self.frames[count], text=text))
                            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                        else:
                            self.frames.append(ttk.Frame(self.Labelframes[4]))
                            self.labels.append(ttk.Label(self.frames[count], text=text))
                            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                
        
        self.pack(fill=BOTH, expand=YES)
        self.scrollframe.pack(fill=BOTH,expand=YES)

        for i in range(len(self.Labelframes)):
            self.Labelframes[i].pack(fill=BOTH, pady=10, padx=20)

        self.update_idletasks()

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)

    def return_values(self):

        m_ivc =self.float_vars(12)
        return  m_ivc, [var.get() for var in self.float_vars]


########### END THIRD SYSTEM #############################################################



############## 4 system MVEM model begin ################################################
  


    
class Notebook_MVEM_model(ttk.Frame):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook, style='cyborg')
    
        self.scrollframe = ScrolledFrame(self)

        self.frames_systems = []

        self.frames = []

        self.Labelframes = [ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Constants"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Initial pressures and temperatures"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Heat capacities"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Coeffcient for pressure drop: delta = k*mdot² [bar*s^2/g^2]"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="CV"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Volume"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Pressure friction coefficients as a function of engine speed p = kf0 + kf1*NE + kf2*NE²"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="air cooler effectiveness: epsilon = kAC0+  kAC1*mdot_a+ kAC2*mdot_a²"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Cylindre"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Different Intertia"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Mass calculations"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Temperature of the air cooler coolant medium [K]")]

        self.labels = []

        self.entrys = []

        self.texts = []

        self.float_vars = []

        self.texts_values= ["pressure air [bar]:", 1,
                            "air Temperature [K]:", 11+273.15,
                            "heat specific capacity ratio for air:",1.401114206,
                            "heat specific capacity ratio for exhaust:", 1.375,
                            "shaft efficiency:",0.99,
                            "combustion efficiency:", 0.99,
                            "turbine efficiency: ", 0.84,
                            "compressor efficiency:", 0.8,
                            
                            "pressure exhaust [bar]:",3.8,
                            "temperature exhuast [K]:",643.15,
                            "temperature inlet [K]:", 298.15,
                            "pressure inlet [bar]:", 2.4,
                            
                            "Air:",1.100,
                            "Exhaust:", 1.006,
                            
                            "k Air cooler:", 3.80*pow(10,-12),
                            "k Air filter:", 0,
                            
                            "cv for inlet:", 0.718,
                            "cv for exhaust:", 0.8,
                            
                            "Inlet receiver [l or dm^3]:", 20.000,
                            "Exhaust receiver [l or dm^3]:",20.000,
                            
                            "kf0:",2.5645,
                            "kf1:",-0.0532,
                            "kf2:",0.0005,
                            
                            "k_ac0:",0.95,
                            "k_ac1:",-5.00*pow(10,-6),
                            "k_ac2:",7*pow(10,-11),

                            "diametre cylindre [m]:",0.92,
                            "number of cylinders :",12,
                            
                            "turbocharger inertia:",600,
                            "engine intertia:",0,
                            "shaft inertia:",0,
                            "propeller inertia:",500000,

                            "coefficient discharge:",1,
                            "Area:",1,
                            "R_a:",287,
                           
                            "T:",32.32+273.15]
        
        
        for i in range(len(self.texts_values)):
            if i%2!=0:
                self.float_vars.append(tk.DoubleVar(value=self.texts_values[i]))
            else:
                self.texts.append(self.texts_values[i])


        for count, text in enumerate(self.texts):
            if count < 8:
                self.frames.append(ttk.Frame(self.Labelframes[0]))
                self.labels.append(ttk.Label(self.frames[count], text=text))
                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
            else:
                if count<12:
                    self.frames.append(ttk.Frame(self.Labelframes[1]))
                    self.labels.append(ttk.Label(self.frames[count], text=text))
                    self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                else:
                    if count<14:
                        self.frames.append(ttk.Frame(self.Labelframes[2]))
                        self.labels.append(ttk.Label(self.frames[count], text=text))
                        self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                    else:
                        if count<16:
                            self.frames.append(ttk.Frame(self.Labelframes[3]))
                            self.labels.append(ttk.Label(self.frames[count], text=text))
                            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                        else:
                            if count<18:
                                self.frames.append(ttk.Frame(self.Labelframes[4]))
                                self.labels.append(ttk.Label(self.frames[count], text=text))
                                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                            else:
                                if count<20:
                                    self.frames.append(ttk.Frame(self.Labelframes[5]))
                                    self.labels.append(ttk.Label(self.frames[count], text=text))
                                    self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                else:
                                    if count<23:
                                        self.frames.append(ttk.Frame(self.Labelframes[6]))
                                        self.labels.append(ttk.Label(self.frames[count], text=text))
                                        self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                    else:
                                        if count<26:
                                            self.frames.append(ttk.Frame(self.Labelframes[7]))
                                            self.labels.append(ttk.Label(self.frames[count], text=text))
                                            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                        else:
                                            if count<28:
                                                self.frames.append(ttk.Frame(self.Labelframes[8]))
                                                self.labels.append(ttk.Label(self.frames[count], text=text))
                                                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                            else:
                                                if count<32:
                                                    self.frames.append(ttk.Frame(self.Labelframes[9]))
                                                    self.labels.append(ttk.Label(self.frames[count], text=text))
                                                    self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                                else:
                                                    if count<35:
                                                        self.frames.append(ttk.Frame(self.Labelframes[10]))
                                                        self.labels.append(ttk.Label(self.frames[count], text=text))
                                                        self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                                    else:
                                                        self.frames.append(ttk.Frame(self.Labelframes[11]))
                                                        self.labels.append(ttk.Label(self.frames[count], text=text))
                                                        self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                    
        
        self.pack(fill=BOTH, expand=YES)
        self.scrollframe.pack(fill=BOTH,expand=YES)

        for i in range(len(self.Labelframes)):
            self.Labelframes[i].pack(fill=BOTH, pady=10, padx=20)

        self.update_idletasks()

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)

    def return_values(self):

        return  [var.get() for var in self.float_vars]
    
############################# END MVEM SYSTEM #############################################################

############################ PARAM Optimisation energy demand BEGIN ################################################

class Notebook_Param_energy(ttk.Frame):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook, style='cyborg')
    
        self.scrollframe = ScrolledFrame(self)

        self.frames_systems = []

        self.frames = []

        self.Labelframes = [ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Constants"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Coefficients for fonction for fuel cons"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Main Engine 1"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Main Engine 2"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Auxiliary Engine 1"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Auxiliary Engine 2")]

        self.labels = []

        self.entrys = []

        self.texts = []

        self.float_vars = []

        self.texts_values= ["Q_ab_max:", 100,
                            "Q_ab_min",3,
                            "Max power",30,
                            "Lower heating value 1:", 100,
                            "Lower heating value 2:", 0.0,
                            "Lower heating value 3:", 0.0,
                            "Lower heating value 4:",0.0 ,
                            "Lower heating value 5:",0.0,
                            "Maximum continuous rating 1:", 63840000,
                            "Maximum continuous rating 2:", 0.0,
                            "Maximum continuous rating 3:", 100,
                            "Maximum continuous rating 4:",0.0 ,
                            "Maximum continuous rating 5:",0.0,

                            "an 1:", 1.0,
                            "an 2:", 1.0,
                            "an 3:", 1.0,
                            "an 4:", 1.0,
                            "an 5:", 1.0,
                            "bn 1:", 1.0,
                            "bn 2:", 1.0,
                            "bn 3:", 1.0,
                            "bn 4:", 1.0,
                            "bn 5:", 1.0,
                            "cn 1:", 1.0,
                            "cn 2:", 1.0,
                            "cn 3:", 1.0,
                            "cn 4:", 1.0,
                            "cn 5:", 1.0,

                            "Maximum power [W]:", 63840000,
                            "Minimum power [W]:", 0.0,
                            "Efficience η1:",0.8 ,
                            "Efficience η2:",0.8,
                            "Efficience ηl:",0.8,

                            "Maximum power [W]:", 1,
                            "Minimum power [W]:", 0.0,
                            "Efficience η1:",0.8 ,
                            "Efficience η2:",0.8,
                            "Efficience ηl:",0.8,

                            "Maximum power [W]:", 1,
                            "Minimum power [W]:", 0.0,
                            "Efficience η1:",0.8 ,
                            "Efficience η2:",0.8,
                            "Efficience ηl:",0.8,

                            "Maximum power [W]:", 1,
                            "Minimum power [W]:", 0.0,
                            "Efficience η1:",0.8 ,
                            "Efficience η2:",0.8,
                            "Efficience ηl:",0.8,]
        
        
        for i in range(len(self.texts_values)):
            if i%2!=0:
                self.float_vars.append(tk.DoubleVar(value=self.texts_values[i]))
            else:
                self.texts.append(self.texts_values[i])


        for count, text in enumerate(self.texts):
            if count < 13:
                # For the ones outside in the first labelframe
                self.frames.append(ttk.Frame(self.Labelframes[0]))
                self.labels.append(ttk.Label(self.frames[count], text=text))
                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
            else:
                if count<28:
                    self.frames.append(ttk.Frame(self.Labelframes[1]))
                    self.labels.append(ttk.Label(self.frames[count], text=text))
                    self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                else:
                    if count<33:
                        self.frames.append(ttk.Frame(self.Labelframes[2]))
                        self.labels.append(ttk.Label(self.frames[count], text=text))
                        self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                    else:
                        if count<38:
                            self.frames.append(ttk.Frame(self.Labelframes[3]))
                            self.labels.append(ttk.Label(self.frames[count], text=text))
                            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                        
                        else:
                            if count<43:
                                self.frames.append(ttk.Frame(self.Labelframes[4]))
                                self.labels.append(ttk.Label(self.frames[count], text=text))
                                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                            else:
                                self.frames.append(ttk.Frame(self.Labelframes[5]))
                                self.labels.append(ttk.Label(self.frames[count], text=text))
                                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                
        
        self.pack(fill=BOTH, expand=YES)
        self.scrollframe.pack(fill=BOTH,expand=YES)

        for i in range(len(self.Labelframes)):
            self.Labelframes[i].pack(fill=BOTH, pady=10, padx=20)

        self.update_idletasks()

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)

    def return_values(self):


        return  [var.get() for var in self.float_vars]