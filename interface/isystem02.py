import tkinter as tk
from tkinter import ttk, filedialog
from ttkbootstrap.constants import *
import pandas as pd


class isystem02(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='cyborg')
        titre_text = "PID Propeller"
        self.titre = ttk.Label(self, text=titre_text, style='main.TLabel')

        self.frames = []

        self.labels = []

        self.entrys = []

        self.texts = ["ship displacement volume [m3]:", "density of sea water [kg/m^3]:", "mass of ship [kg]:",
                      "add virtual mass [kg]:", "diameter [m]:", "thrust:", "dVs_dt = c_1 * V_s^2 [kg/m]:",
                      "# the number of propeller blades, zp:", "the disk area coefficient, AE/Ao:",
                      "the pitch to diameter ratio, p/Dp:", "ship wake fraction, w, which is considered constant "
                      "taking values in the range from 0.20:"]
        
        self.float_vars = [tk.DoubleVar() for _ in self.texts]

        for count, text in enumerate(self.texts):
            self.frames.append(ttk.Frame(self))
            self.labels.append(ttk.Label(self.frames[count], text=text))
            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count]))

        # Adiciona um traço divisor
        self.divider1 = ttk.Separator(self, orient='horizontal')
        self.divider2 = ttk.Separator(self, orient='horizontal')

        # Campo de arquivo para os coeficientes de Wageningen
        self.wageningen_frame = ttk.Frame(self)
        self.wageningen_label = ttk.Label(self.wageningen_frame, text="Wageningen coefficients imports:")
        self.wageningen_file_entry = ttk.Entry(self.wageningen_frame, state='readonly')
        self.browse_button = ttk.Button(self.wageningen_frame, text="Browse", command=self.browse_wageningen_file)

    def show(self):
        self.pack(fill=BOTH, expand=YES)
        self.titre.pack(fill=X, pady=10, padx=20)
        self.update_idletasks()

        for i in range(len(self.labels)):
            self.frames[i].pack(side=TOP, anchor='w')
            self.labels[i].pack(side=LEFT, padx=5, pady=5)
            self.entrys[i].pack(side=LEFT, padx=5, pady=5)
            if i==1:
                self.divider1.pack(fill=X, pady=5)


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
        file_path = filedialog.askopenfilename(filetypes=[("Excel files", "*.xlsx;*.xls")])
        self.wageningen_file_entry.config(state='normal')
        self.wageningen_file_entry.delete(0, tk.END)
        self.wageningen_file_entry.insert(0, file_path)
        self.wageningen_file_entry.config(state='readonly')
        # Aqui você pode adicionar a lógica para lidar com o arquivo Excel carregado
        # Exemplo: leitura dos coeficientes do arquivo e atualização das variáveis de controle
        df = pd.read_excel(file_path)
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
        
        print(f"Arquivo Excel selecionado: {file_path}")

