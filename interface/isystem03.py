import tkinter as tk
import ttkbootstrap as ttk

from ttkbootstrap.constants import *


class isystem03(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='cyborg')
        titre_text = "Powermanagement"
        self.titre = ttk.Label(self, text=titre_text, style='main.TLabel')
        
        self.float_vars = [tk.DoubleVar() for _ in range(7)]

        self.frames = []

        self.labels = []

        self.entrys = []

        self.texts = ["Input 1:","Input 2:","Input 3:","Input 4:",
                      "Input 5:","Input 6:"]
        
        for count,text in enumerate(self.texts):

            self.frames.append(ttk.Frame(self))

            self.labels.append(ttk.Label(self.frames[count], text= text))

            self.entrys.append(ttk.Entry(self.frames[count], textvariable=self.float_vars[count])) 

            


    def show(self):
        self.pack(fill=BOTH, expand=YES)
        self.titre.pack(fill=X, pady=10, padx=20)
        self.update_idletasks() 
        
        for i in range(len(self.labels)):
                self.frames[i].pack()
                self.labels[i].pack(side=LEFT, padx=5, pady=5)
                self.entrys[i].pack(side=RIGHT, padx=5, pady=5)
            
            
        

    def unshow(self):
        self.pack_forget()
    


