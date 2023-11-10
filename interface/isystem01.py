import tkinter as tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *


class isystem01(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='secondary')
        self.pack(side=LEFT,fill=Y)

        system01 = ttk.Button(self, text='Perra',width=15)
        
        

