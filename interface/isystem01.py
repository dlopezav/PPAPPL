import tkinter as tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *




class isystem01(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='cyborg')
        titre_text = "System 1"

        self.titre = ttk.Label(self, text=titre_text, style='main.TLabel')
        
        self.float_vars = [tk.DoubleVar() for _ in range(7)]
        
        self.labels = [ttk.Label(self, text="Input 1:"),ttk.Label(self, text="Input 2:"),
                       ttk.Label(self, text="Input 3:"),ttk.Label(self, text="Input 4:"),
                       ttk.Label(self, text="Input 5:"),ttk.Label(self, text="Input 6:"),]

        self.entrys = [ttk.Entry(self, textvariable=self.float_vars[0]),ttk.Entry(self, textvariable=self.float_vars[1]),
                       ttk.Entry(self, textvariable=self.float_vars[2]),ttk.Entry(self, textvariable=self.float_vars[3]),
                       ttk.Entry(self, textvariable=self.float_vars[4]),ttk.Entry(self, textvariable=self.float_vars[5])]

            


    def show(self):
        self.pack(fill=BOTH, expand=YES)
        self.titre.pack(fill=X, pady=10, padx=20)
        self.update_idletasks() 
        
        for i in range(len(self.labels)):
            self.labels[i].pack(side=TOP, padx=5, pady=5)             
            self.entrys[i].pack(side=TOP, padx=5, pady=5)

            
        

    def unshow(self):
        self.pack_forget()
    


