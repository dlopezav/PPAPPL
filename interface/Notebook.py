

###################################### Frame classes to add in main Notebook  #########################################
# Libraries needed
import tkinter as tk
from tkinter import ttk, filedialog
from ttkbootstrap.constants import *
from ttkbootstrap.scrolled import ScrolledFrame
import pandas as pd
import math
import openpyxl


######################################## MOTHER CLASS PAGE #######################################################
# First mother class Notebook_page that represents a general tab, initializing the variables 
# This page contains the main methods

class Notebook_page(ttk.Frame):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook, style='cyborg')

        self.frames_systems = []

        self.Labelframes = []

        self.frames = []
        
        self.scrollframe = ScrolledFrame(self)

        self.file_path= ""

        self.labels = []

        self.entrys = []

        self.texts = []

        self.reset = []

        self.float_vars = []

        self.texts_values= []

        self.files = [""]

        self.file_entry = []
    
    # Method to browse an excel file from the diretory
    def browse_file(self, number):
        self.file_path = filedialog.askopenfilename(initialdir="../",
                                                    filetypes=[("Excel files", ".xlsx .xls")],
                                                    title="Choose a file.")

        self.file_entry[number].config(state='normal')
        self.file_entry[number].delete(0, tk.END)
        self.file_entry[number].insert(0, self.file_path)
        self.file_entry[number].config(state='readonly')
        self.files[number] =  self.file_path
    
    # Method to choose a diretory to save plots and steady states sheets
    def browse_dir(self, number):
        self.dir_path = filedialog.askdirectory(initialdir="../",
                                                title="Choose a directory.")
        self.file_entry[number].config(state='normal')
        self.file_entry[number].delete(0, tk.END)
        self.file_entry[number].insert(0, self.dir_path)
        self.file_entry[number].config(state='readonly')
        self.files[number] =  self.dir_path
    
    def browse_csv(self, number):
        self.file_path = filedialog.askopenfilename(initialdir="../",
                                                    filetypes=[("CSV files", "*.csv")],
                                                    title="Choose a file.")

        self.file_entry[number].config(state='normal')
        self.file_entry[number].delete(0, tk.END)
        self.file_entry[number].insert(0, self.file_path)
        self.file_entry[number].config(state='readonly')
        self.files[number] =  self.file_path

    # Method to return all the values in the system
    def return_values(self):

        if len(self.files)>1:
            return self.files
        else:
            return  [var.get() for var in self.float_vars]
    
    # Method to reset values to their default values
    def setDefaultVal(self, i):
        self.entrys[i].delete(0, tk.END)
        self.entrys[i].insert(0, self.default[i*2+1])

    def charge_paths(self):
        if len(self.file_entry) != 0:
            for number in range(len(self.files)):
                self.file_entry[number].config(state='normal')
                self.file_entry[number].delete(0, tk.END)
                self.file_entry[number].insert(0, self.files[number])
                self.file_entry[number].config(state='readonly')
        



################################ MODEL PARAMETERS BEGIN ###################################################
'''
They following classes are subclasses of the notebook tab
They have all the same structure
Some of them has files as an input

'''

# Define a class Notebook_Model_param that inherits from Notebook_page
class Notebook_Model_param(Notebook_page):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook)

        # Initialize label frames for the notebook
        self.Labelframes = [ttk.LabelFrame(self, bootstyle='info', text="MODEL")]

        # Define default values for various parameters
        self.texts_values = ["time_t:", 0,  # time_t 0 
                            "total_time:", 0.005 ,  # total_time 1
                            "deltat2:",0.005,  # delta2 2
                            "initial value of turbocharger speed [unit?]:", 7522,  # NTC_ini 3
                            "initial pressure ratio compressor:",2.4,  # prc 4
                            "initial pressure ratio turbine:", 3.8,  # prt 5
                            "Initial_speed: ", 12.15,  # Vs_ini 6
                            "model time of the simulation:", 0.035,  # total_time 7
                            ]

        # Store default values for resetting
        self.default = self.texts_values.copy()

        # Initialize lists for storing labels, entry widgets, and reset buttons
        for i in range(len(self.texts_values)):
            if i % 2 != 0:
                self.float_vars.append(tk.DoubleVar(value=self.texts_values[i]))
            else:
                self.texts.append(self.texts_values[i])

        # Create frames, labels, entry widgets, and reset buttons
        for count, text in enumerate(self.texts):
            self.frames.append(ttk.Frame(self.Labelframes[0]))
            self.labels.append(ttk.Label(self.frames[count], text=text))
            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
            self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count: self.setDefaultVal(count), style='danger-link'))

        # Pack widgets and label frames
        self.pack(fill=BOTH, expand=YES)

        for i in range(len(self.Labelframes)):
            self.Labelframes[i].pack(fill=BOTH, pady=10, padx=20)

        self.update_idletasks()

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)
            self.reset[i].pack(side=RIGHT, padx=5, pady=5)

        # Calculate and append an additional value to float_vars (N_cycle)
        self.float_vars.append(tk.DoubleVar(value=int(self.float_vars[1].get() // self.float_vars[2].get() + 1)))  # N_cycle 8

    
################ First system propeller ################################

class Notebook_propeller(Notebook_page):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook)

        self.Labelframes = [ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Constants"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="To be determined")]

        self.texts_values =["Ship displacement volume [m3]:", 112404,# sh_displ 9
                            "Density of sea water [kg/m^3]:", 1026,# rho_sw 10
                            "Mass of ship [kg]:", 260000000,# mass_of_ship 11 
                            "Diameter [m]:", 10,# D_p 12 
                            "Thrust:", 0.2,# thrust 13
                            "dVs_dt = c_1 * V_s^2 [kg/m]:", 26250,# c_1 14
                            "Number of propeller blades, zp:", 6,# z_p 15
                            "Disk area coefficient, AE/Ao:", 0.82,# A_e_A_o 16
                            "Pitch to diameter ratio, p/Dp:", 0.93846154,# p_D_p 17
                            "Ship wake fraction, w, which is considered constant \n taking values in the range from 0.20:", 0.3 #omega 18
                            ]
        
        self.default = self.texts_values.copy()

        


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
                self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))
            else:
                self.frames.append(ttk.Frame(self.Labelframes[1]))
                self.labels.append(ttk.Label(self.frames[count], text=text))
                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda  count=count: self.setDefaultVal(count), style='danger-link'))

        
        # Wageningen
        self.file_frame = ttk.Frame(self.scrollframe)
        self.file_label = ttk.Label(self.file_frame, text="Wageningen coefficients imports:")
        self.file_entry.append(ttk.Entry(self.file_frame, state='readonly'))
        self.file_button = ttk.Button(self.file_frame, text="Browse", command=lambda : self.browse_file(0), style='info')
        
        self.pack(fill=BOTH, expand=YES)
        self.scrollframe.pack(fill=BOTH, expand=YES)   
        self.Labelframes[0].pack(fill=BOTH, pady=10, padx=20)
        self.Labelframes[1].pack(fill=BOTH, pady=10, padx=20)
        self.update_idletasks()

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)
            self.reset[i].pack(side=RIGHT, padx=5, pady=5)

        
        # Position wageningen frame and label
        self.file_frame.pack(side=TOP, pady=5)
        self.file_label.pack(side=TOP, padx=5)
        self.file_button.pack(side=LEFT, padx=5)
        self.file_entry[0].pack(side=LEFT, padx=5)
        self.count=FALSE

        if (self.file_path!="" and count==FALSE):
            df = pd.read_excel(self.file_path)

            self.float_vars.append(tk.DoubleVar(value=float(df['n1'].iloc[0])))#n1_nump 19
            self.float_vars.append(tk.DoubleVar(value=float(df['ct'].iloc[0])))#ct_nump 20
            self.float_vars.append(tk.DoubleVar(value=float(df['s'].iloc[0])))#s1_nump 21
            self.float_vars.append(tk.DoubleVar(value=float(df['t'].iloc[0])))#t1_nump 22 
            self.float_vars.append(tk.DoubleVar(value=float(df['u'].iloc[0])))#u1_nump 23
            self.float_vars.append(tk.DoubleVar(value=float(df['v'].iloc[0])))#v1_nump 24

            self.float_vars.append(tk.DoubleVar(value=float(df['n2'].iloc[0])))#n2_nump 25
            self.float_vars.append(tk.DoubleVar(value=float(df['cq'].iloc[0])))#cq_nump 26
            self.float_vars.append(tk.DoubleVar(value=float(df['s2'].iloc[0])))#s2_nump 27 
            self.float_vars.append(tk.DoubleVar(value=float(df['t2'].iloc[0])))#t2_nump 28
            self.float_vars.append(tk.DoubleVar(value=float(df['u2'].iloc[0])))#u2_nump 29
            self.float_vars.append(tk.DoubleVar(value=float(df['v2'].iloc[0])))#v2_nump 30

            self.float_vars.append(tk.DoubleVar(value=float(self.float_vars[0].get() * self.float_vars[1].get())))#m_hydro 31
            self.count=TRUE

    

################ END FIRST SYSTEM ######################################################

################## Second system PID controller #########################################

class Notebook_PIDcontroller(Notebook_page):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook)

        self.Labelframes = [ttk.LabelFrame(self, bootstyle='info',text="Constants"),
                            ttk.LabelFrame(self, bootstyle='info',text="To be determined")]

        self.texts_values =["Number of engine cylinders:", 12, #zc 32
                            "Constant for PID controller	kp:", 0.003, #kp 33
                            "Constant for PID controller	kd:", 0, #kd  34
                            "Constant for PID controller	ki:", 0.0015, #ki 35
                            "Original injection time (at the beginning of the simulation):",0.01427586 #xr_o 36
                            ]
        
        

        self.default = self.texts_values.copy()
        
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
                self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

            else:
                self.frames.append(ttk.Frame(self.Labelframes[1]))
                self.labels.append(ttk.Label(self.frames[count], text=text))
                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))


        
        # Steadystate
        self.file_frame = ttk.Frame(self)
        self.file_label = ttk.Label(self.file_frame, text="Steadystate imports:")
        self.file_entry.append(ttk.Entry(self.file_frame, state='readonly'))
        self.file_button = ttk.Button(self.file_frame, text="Browse", command=lambda :self.browse_file(0), style='info')
        
        self.pack(fill=BOTH, expand=YES)   
        self.Labelframes[0].pack(fill=BOTH, pady=10, padx=20)
        self.Labelframes[1].pack(fill=BOTH, pady=10, padx=20)
        self.update_idletasks()
    

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)
            self.reset[i].pack(side=RIGHT, padx=5, pady=5)
                
        

        # Position steadystate frame and label
        self.file_frame.pack(side=TOP, pady=5)
        self.file_label.pack(side=TOP, padx=5)
        self.file_button.pack(side=LEFT, padx=5)
        self.file_entry[0].pack(side=LEFT, padx=5)
        self.count = FALSE

  
        if (self.file_path!="" and self.count == FALSE):

            steadystate = pd.read_excel(self.file_path)
            me_ne_value = steadystate['me_ne'].iloc[0]
            calculated_value = 2 * math.pi / 60 * me_ne_value
            self.float_vars.append(tk.DoubleVar(value=float(calculated_value))) #Nord 37
            self.count = TRUE


########### END SECOND SYSTEM #############################################################



############## Third system 0d cycle begin ################################################

class Notebook_0dCycle(Notebook_page):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook)


        self.Labelframes = [ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Constants"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="To be confirmed"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Wiebe parameters"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Parameters combustion p167"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Engine reference point")]


        self.texts_values= ["Delta t:", 0.0001,#deltat 38
                            "Number of revolution per cycle:", 1,#rev_cy 39
                            "Clearance volume [cm3]:",164671,#Vc 40
                            "Lower heating value :", 42707,#lhv 41
                            "wall surface [units ?]:",1.102343986,#Sw 42
                            "Stroke:", 3.468,#stroke 43
                            "l: ", 3468,#l 44
                            "ar:", 1734,#ar 45
                            "Perfect gas constant:",8.314462,#R 46

                            "Tw [K]:",300+273.15,#Tw 47
                            "Pressure at inlet valve closing [bar]:", 4.1,#p_ivc 48
                            "Temperature at the inlet valve closing [K]:", 307,#T_ivc 49
                            "Initial masse in the cylinder [g]:",7462,#ma_ini 50
                            "Compression ratio:", 15,#rc 51
                            "reference temperature [K]:", 298.15,#Tref 52
                            "Stochiometric coefficient:", 14.8,#Cstoich 53
                            "R / molar mass of exhaust", 0.1341,#re 54
                            "R / mass of air", 0.287,#ra 55
                            "R / mass of fuel gas", 0.0489086,#rf 56
                            
                            "acomb:",0.42,#acomb 57
                            "bcomb:",0.49,#bcomb 58
                            "am:",0.59,#am 59
                            "bm:",0.21,#bm 60
                            "cm:",-0.96,#cm 61
                            "delta_m:",0.27,#delta_m 62

                            "a_id:",0.39,#a_id 63
                            "b_id:",0.105,#b_id 64
                            "c_id:",3.12,#c_id 65
                            "mWiebe:",0.1,#mWiebe 66

                            "Air fuel ratio for reference point:",2.435696566131349,#lamb_ref 67
                            "ignition delay for reference point",1,#phi_id_ref 68
                            "Mass at inlet valve closing for reference point [g]:",7462,#m_ivc_ref 69
                            "Engine speed for reference point [rpm]:",80,#N_re 70
                            "Wiebe curve form coefficient for reference point:",0.1,#mWieberef 71
                            "Combustion duration for reference point:",0.3667076110839844*60/(2*math.pi),#delta_phi_comb_ref 72
                            "SMFR:",14500,#SMFR 73
                            "topen:",0.0005,#topen 74
                            "tclose",0.0005#tclose 75
                            ]
        
        self.default = self.texts_values.copy()
        
        
        
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
                self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

            else:
                if count<19:
                    self.frames.append(ttk.Frame(self.Labelframes[1]))
                    self.labels.append(ttk.Label(self.frames[count], text=text))
                    self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                    self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                else:
                    if count<25:
                        self.frames.append(ttk.Frame(self.Labelframes[2]))
                        self.labels.append(ttk.Label(self.frames[count], text=text))
                        self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                        self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                    else:
                        if count<29:
                            self.frames.append(ttk.Frame(self.Labelframes[3]))
                            self.labels.append(ttk.Label(self.frames[count], text=text))
                            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                            self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                        else:
                            self.frames.append(ttk.Frame(self.Labelframes[4]))
                            self.labels.append(ttk.Label(self.frames[count], text=text))
                            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                            self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                
        

        self.pack(fill=BOTH, expand=YES)
        self.scrollframe.pack(fill=BOTH,expand=YES)

        for i in range(len(self.Labelframes)):
            self.Labelframes[i].pack(fill=BOTH, pady=10, padx=20)

        self.update_idletasks()
        


        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)
            self.reset[i].pack(side=RIGHT, padx=5, pady=5)


########### END THIRD SYSTEM #############################################################



############## 4 system MVEM model begin ################################################
  


    
class Notebook_MVEM_model(Notebook_page):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook)

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


        self.texts_values= ["pressure air [bar]:", 1,#p_a 76
                            "air Temperature [K]:", 11+273.15,#T_a 77
                            "heat specific capacity ratio for air:",1.401114206,#gamma_a 78
                            "heat specific capacity ratio for exhaust:", 1.375,#gamma_e 79 
                            "shaft efficiency:",0.99,#eta_sh 80 
                            "combustion efficiency:", 0.99,#eta_comb 81 
                            "turbine efficiency: ", 0.84,#eta_t 82
                            "compressor efficiency:", 0.8,#eta_c 83
                            
                            "pressure exhaust [bar]:",3.8,#p_ER_ini 84
                            "temperature exhuast [K]:",643.15,#T_ER_ini 85
                            "temperature inlet [K]:", 298.15,#T_IR_ini 86
                            "pressure inlet [bar]:", 2.4,#p_IR_ini 87
                            
                            "Air:",1.100,#cpe 88
                            "Exhaust:", 1.006,#cpa 89
                            
                            "k Air cooler:", 3.80*pow(10,-12),#kpac 90
                            "k Air filter:", 0,#kpaf 91
                            
                            "cv for inlet:", 0.718,#cvir 92
                            "cv for exhaust:", 0.8,#cver 93
                            
                            "Inlet receiver [l or dm^3]:", 20.000,#V_ir 94
                            "Exhaust receiver [l or dm^3]:",20.000,#V_er 95
                            
                            "kf0:",2.5645,#kf0 96
                            "kf1:",-0.0532,#kf1 97 
                            "kf2:",0.0005,#kf2 98
                            
                            "k_ac0:",0.95,#k_ac0 99
                            "k_ac1:",-5.00*pow(10,-6),#k_ac1 100
                            "k_ac2:",7*pow(10,-11),#k_ac2 101

                            "diametre cylindre [m]:",0.92,#bore 102
                            "number of cylinders :",12,#zc 103
                            
                            "turbocharger inertia:",600,#ITC 104
                            "engine intertia:",0,#I_e 105
                            "shaft inertia:",0,#I_sh 106
                            "propeller inertia:",500000,#I_p 107

                            "coefficient discharge:",1,#cd 108
                            "Area:",1,#A_eq 109
                            "R_a:",287,#R_a 110
                           
                            "T:",32.32+273.15#T_w 111
                            ]
        
        self.default = self.texts_values.copy()
        
        
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
                self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

            else:
                if count<12:
                    self.frames.append(ttk.Frame(self.Labelframes[1]))
                    self.labels.append(ttk.Label(self.frames[count], text=text))
                    self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                    self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                else:
                    if count<14:
                        self.frames.append(ttk.Frame(self.Labelframes[2]))
                        self.labels.append(ttk.Label(self.frames[count], text=text))
                        self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                        self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                    else:
                        if count<16:
                            self.frames.append(ttk.Frame(self.Labelframes[3]))
                            self.labels.append(ttk.Label(self.frames[count], text=text))
                            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                            self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                        else:
                            if count<18:
                                self.frames.append(ttk.Frame(self.Labelframes[4]))
                                self.labels.append(ttk.Label(self.frames[count], text=text))
                                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                            else:
                                if count<20:
                                    self.frames.append(ttk.Frame(self.Labelframes[5]))
                                    self.labels.append(ttk.Label(self.frames[count], text=text))
                                    self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                    self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                                else:
                                    if count<23:
                                        self.frames.append(ttk.Frame(self.Labelframes[6]))
                                        self.labels.append(ttk.Label(self.frames[count], text=text))
                                        self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                        self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                                    else:
                                        if count<26:
                                            self.frames.append(ttk.Frame(self.Labelframes[7]))
                                            self.labels.append(ttk.Label(self.frames[count], text=text))
                                            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                            self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                                        else:
                                            if count<28:
                                                self.frames.append(ttk.Frame(self.Labelframes[8]))
                                                self.labels.append(ttk.Label(self.frames[count], text=text))
                                                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                                self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                                            else:
                                                if count<32:
                                                    self.frames.append(ttk.Frame(self.Labelframes[9]))
                                                    self.labels.append(ttk.Label(self.frames[count], text=text))
                                                    self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                                    self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                                                else:
                                                    if count<35:
                                                        self.frames.append(ttk.Frame(self.Labelframes[10]))
                                                        self.labels.append(ttk.Label(self.frames[count], text=text))
                                                        self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                                        self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                                                    else:
                                                        self.frames.append(ttk.Frame(self.Labelframes[11]))
                                                        self.labels.append(ttk.Label(self.frames[count], text=text))
                                                        self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                                        self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                                    
        

        self.pack(fill=BOTH, expand=YES)
        self.scrollframe.pack(fill=BOTH,expand=YES)

        for i in range(len(self.Labelframes)):
            self.Labelframes[i].pack(fill=BOTH, pady=10, padx=20)

        self.update_idletasks()
      

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)
            self.reset[i].pack(side=RIGHT, padx=5, pady=5)
        
        

    
############################# END MVEM SYSTEM #############################################################

############################ PARAM Optimisation energy demand BEGIN ################################################

class Notebook_Param_energy(Notebook_page):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook)

        self.Labelframes = [ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Constants"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Coefficients for fonction for fuel cons"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Main Engine 1"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Main Engine 2"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Auxiliary Engine 1"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Auxiliary Engine 2")]

        self.texts_values= ["Q_ab_max:", 100,#Q_ab_max 112
                            "Q_ab_min",3,#Q_ab_min 113
                            "Max power",30,#Max_power 114
                            "Lower heating value 1:", 100,#LHV1 115
                            "Lower heating value 2:", 0.0,#LHV2 116
                            "Lower heating value 3:", 0.0,#LHV3 117
                            "Lower heating value 4:",0.0 ,#LHV4 118
                            "Lower heating value 5:",0.0,#LHV5 119
                            "Maximum continuous rating 1:", 63840000,#MCR1 120
                            "Maximum continuous rating 2:", 0.0,#MCR2 121
                            "Maximum continuous rating 3:", 100,#MCR3 122
                            "Maximum continuous rating 4:",0.0 ,#MCR4 123
                            "Maximum continuous rating 5:",0.0,#MCR5 124

                            "an 1:", 1.0,#an 125
                            "an 2:", 1.0,#an2 126
                            "an 3:", 1.0,#an3 127
                            "an 4:", 1.0,#an4 128
                            "an 5:", 1.0,#an5 129
                            "bn 1:", 1.0,#bn 130
                            "bn 2:", 1.0,#bn2 131
                            "bn 3:", 1.0,#bn3 132
                            "bn 4:", 1.0,#bn4 133
                            "bn 5:", 1.0,#bn5 134
                            "cn 1:", 1.0,#cn 135
                            "cn 2:", 1.0,#cn2 136
                            "cn 3:", 1.0,#cn3 137
                            "cn 4:", 1.0,#cn4 138
                            "cn 5:", 1.0,#cn5 139
 
                            "Maximum power [W]:", 63840000,#PME1_max 140
                            "Minimum power [W]:", 0.0,#PME1_min 141 
                            "Efficience η1:",0.8 ,#nu1ME1 142
                            "Efficience η2:",0.8,#nu2ME1 143
                            "Efficience ηl:",0.8,#nuelME1 144
 
                            "Maximum power [W]:", 1,#PME2_max 145
                            "Minimum power [W]:", 0.0,#PME2_min 146
                            "Efficience η1:",0.8 ,#nu1ME2 147
                            "Efficience η2:",0.8,#nu2ME2 148
                            "Efficience ηl:",0.8,#nuelME2 149

                            "Maximum power [W]:", 1,#PAE1_max 150
                            "Minimum power [W]:", 0.0,#PAE1_min 151
                            "Efficience η1:",0.8 ,#nu1AE1 152
                            "Efficience η2:",0.8,#nu2AE1 153
                            "Efficience ηl:",0.8,#nuelAE1 154
 
                            "Maximum power [W]:", 1,#PAE2_max 155
                            "Minimum power [W]:", 0.0,#PAE2_min 156
                            "Efficience η1:",0.8 ,#nu1AE2 157
                            "Efficience η2:",0.8,#nu2AE2 158
                            "Efficience ηl:",0.8,#nuelAE2 159
                            ]
        
        self.default = self.texts_values.copy()

        
        
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
                self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

            else:
                if count<28:
                    self.frames.append(ttk.Frame(self.Labelframes[1]))
                    self.labels.append(ttk.Label(self.frames[count], text=text))
                    self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                    self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                else:
                    if count<33:
                        self.frames.append(ttk.Frame(self.Labelframes[2]))
                        self.labels.append(ttk.Label(self.frames[count], text=text))
                        self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                        self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                    else:
                        if count<38:
                            self.frames.append(ttk.Frame(self.Labelframes[3]))
                            self.labels.append(ttk.Label(self.frames[count], text=text))
                            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                            self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                        
                        else:
                            if count<43:
                                self.frames.append(ttk.Frame(self.Labelframes[4]))
                                self.labels.append(ttk.Label(self.frames[count], text=text))
                                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                            else:
                                self.frames.append(ttk.Frame(self.Labelframes[5]))
                                self.labels.append(ttk.Label(self.frames[count], text=text))
                                self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))
                                self.reset.append(ttk.Button(self.frames[count], text="Reset value", command=lambda count=count : self.setDefaultVal(count), style='danger-link'))

                
       

        self.pack(fill=BOTH, expand=YES)
        self.scrollframe.pack(fill=BOTH,expand=YES)

        for i in range(len(self.Labelframes)):
            self.Labelframes[i].pack(fill=BOTH, pady=10, padx=20)

        self.update_idletasks()
      

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)
            self.reset[i].pack(side=RIGHT, padx=5, pady=5)
        
        

###################################### EXCEL SHEETS #####################################################

class Notebook_excell_sheets(Notebook_page):
    def __init__(self, tab_notebook):
        super().__init__(tab_notebook)


        self.Labelframes = [ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Steady states engines 25,30...70"),
                            ttk.LabelFrame(self.scrollframe, bootstyle='info',text="Others (excel_files)")]

        self.texts_values =["Choose folder for steady_states_engines_xx sheets",
                            "Choose directory to save plots:", 
                            

                            "Results: ",  
                            "TC_maps2:", 
                            "Wageningen_reynolds:",
                            "KT:", 
                            "KQ:", 
                            "steady_states:",
                            "Directory for plotting results over selected cycles for steady state",
                            "Pressure real vs cantera",
                            "trial_DF25"
                            ]
        
        self.files = self.files*len(self.texts_values)
        
        self.buttons = []
        
        for count, text in enumerate(self.texts_values):
            if count<2:
                self.frames.append(ttk.Frame(self.Labelframes[0]))
                self.labels.append(ttk.Label(self.frames[count], text=text))
                self.file_entry.append(ttk.Entry(self.frames[count], state='readonly'))
                self.buttons.append(ttk.Button(self.frames[count], text="Browse", command=lambda c=count: self.browse_dir(c), style='info'))
            else:
                if count<8:
                    self.frames.append(ttk.Frame(self.Labelframes[1]))
                    self.labels.append(ttk.Label(self.frames[count], text=text))
                    self.file_entry.append(ttk.Entry(self.frames[count], state='readonly'))
                    self.buttons.append(ttk.Button(self.frames[count], text="Browse", command=lambda c=count: self.browse_file(c), style='info'))
                else:
                    if count==8:
                        self.frames.append(ttk.Frame(self.Labelframes[1]))
                        self.labels.append(ttk.Label(self.frames[count], text=text))
                        self.file_entry.append(ttk.Entry(self.frames[count], state='readonly'))
                        self.buttons.append(ttk.Button(self.frames[count], text="Browse", command=lambda c=count: self.browse_dir(c), style='info'))
                    else:
                        self.frames.append(ttk.Frame(self.Labelframes[1]))
                        self.labels.append(ttk.Label(self.frames[count], text=text))
                        self.file_entry.append(ttk.Entry(self.frames[count], state='readonly'))
                        self.buttons.append(ttk.Button(self.frames[count], text="Browse", command=lambda c=count: self.browse_csv(c), style='info'))

                

        
        self.pack(fill=BOTH, expand=YES)  
        self.scrollframe.pack(fill=BOTH, expand=YES) 
        self.Labelframes[0].pack(fill=BOTH, pady=10, padx=20)
        self.Labelframes[1].pack(fill=BOTH, pady=10, padx=20)
        self.update_idletasks()
        

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.file_entry[i].pack(side=LEFT, padx=5, pady=5)
            self.buttons[i].pack(side=LEFT, padx=5)


        
