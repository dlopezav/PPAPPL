import tkinter as tk
from tkinter import ttk, filedialog
from ttkbootstrap.constants import *
from ttkbootstrap.scrolled import ScrolledFrame
import pandas as pd
import math

class main_page(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='cyborg')

        titre_text = "MAIN"

        self.titre = ttk.Label(self, text=titre_text, style='main.TLabel')

        self.values = []

        self.notebook = ttk.Notebook(self, bootstyle="info")

        self.frames_notebook = [Notebook_propeller(self.notebook),Notebook_PIDcontroller(self.notebook),
                                Notebook_0dCycle(self.notebook)]

        self.notebook.add(self.frames_notebook[0],text="PID propeller")
        self.notebook.add(self.frames_notebook[1],text="PID controller")
        self.notebook.add(self.frames_notebook[2],text="od Cycle")

    def show(self):
        self.pack(fill=BOTH, expand=YES)
        self.titre.pack(fill=X, pady=10, padx=20)
        self.notebook.pack(fill=BOTH,expand=YES)
        self.update_idletasks()


    def unshow(self):
        self.pack_forget()

    def return_all_values(self):
        for i in range(len(self.frames_notebook)):
            self.values.append(self.frames_notebook[i].return_values())
        
        return self.values
        

###################################### Frame classes to add in the Notebook  #########################################


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

        self.texts = ["Ship displacement volume [m3]:", "Density of sea water [kg/m^3]:", "Mass of ship [kg]:",
                      "Add virtual mass [kg]:", "Diameter [m]:", "Thrust:", "dVs_dt = c_1 * V_s^2 [kg/m]:",
                      "Number of propeller blades, zp:", "Disk area coefficient, AE/Ao:",
                      "Pitch to diameter ratio, p/Dp:", "Ship wake fraction, w, which is considered constant "
                      "Taking values in the range from 0.20:"]
        
        self.float_vars = [tk.DoubleVar() for _ in self.texts]

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

        return n1_nump, ct_nump, s1_nump, t1_nump, u1_nump, v1_nump, n2_nump, cq_nump, s2_nump, t2_nump, u2_nump, v2_nump, \
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

        self.texts = ["Number of engine cylinders:", 
                      "Constant for PID controller	kp:", 
                      "Constant for PID controller	kd:",
                      "Constant for PID controller	ki:", 
                      "Original injeciton time (at the beginning of the simulation):"]
        
        self.float_vars = [tk.DoubleVar(),tk.DoubleVar(value=0.003),
                           tk.DoubleVar(),tk.DoubleVar(value=0.0015),tk.DoubleVar(value=0.01427586)]

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


        return  [var.get() for var in self.float_vars]