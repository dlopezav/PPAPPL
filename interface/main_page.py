import tkinter as tk
from tkinter import ttk, filedialog
from ttkbootstrap.constants import *
from ttkbootstrap.scrolled import ScrolledFrame
import pandas as pd

class main_page(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='cyborg')

        titre_text = "MAIN"

        self.titre = ttk.Label(self, text=titre_text, style='main.TLabel')

        self.notebook = ttk.Notebook(self, bootstyle="cyborg")

        self.frame_to_scroll = ttk.Frame(self.notebook)

        self.notebook.add(self.frame_to_scroll,text="propeller")

        self.scrollframe = ScrolledFrame(self.frame_to_scroll)

        self.notebook.add(self, text="propeller")
    
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

        # Separator
        self.divider1 = ttk.Separator(self, orient='horizontal')
        self.divider2 = ttk.Separator(self, orient='horizontal')

        # Wageningen
        self.wageningen_frame = ttk.Frame(self.scrollframe)
        self.wageningen_label = ttk.Label(self.wageningen_frame, text="Wageningen coefficients imports:")
        self.wageningen_file_entry = ttk.Entry(self.wageningen_frame, state='readonly')
        self.browse_button = ttk.Button(self.wageningen_frame, text="Browse", command=self.browse_wageningen_file)
        

    def show(self):
        self.pack(fill=BOTH, expand=YES)
        self.titre.pack(fill=X, pady=10, padx=20)
        self.notebook.pack(fill=BOTH,expand=YES)
        self.scrollframe.pack(fill=BOTH, expand=YES)   
        self.Labelframes[0].pack(fill=BOTH, pady=10, padx=20)
        self.Labelframes[1].pack(fill=BOTH, pady=10, padx=20)
        self.update_idletasks()

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)
                


        # Adiciona o traço divisor
        self.divider2.pack(fill=X, pady=5)

        # Posicionamento do campo de arquivo próximo ao parâmetro de Wageningen
        self.wageningen_frame.pack(side=TOP, pady=5)
        self.wageningen_label.pack(side=TOP, padx=5)
        self.browse_button.pack(side=LEFT, padx=5)
        self.wageningen_file_entry.pack(side=LEFT, padx=5)

        

    def unshow(self):
        self.pack_forget()

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
