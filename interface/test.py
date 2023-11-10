import tkinter as tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from ttkbootstrap import Style


class MainProgram(ttk.Frame):
    def __init__(self, master_window):
        super().__init__(master_window)
        self.pack(fill=BOTH, expand=YES)

        NavBar(self)

        self.frames = {}
        
        for F
        for F in (StartPage, Page1, Page2):

			frame = F(container, self)

			# initializing frame of that object from
			# startpage, page1, page2 respectively with 
			# for loop
			self.frames[F] = frame 

			frame.grid(row = 0, column = 0, sticky ="nsew")

		self.show_frame(StartPage)

	# to display the current frame passed as
	# parameter
	def show_frame(self, cont):
		frame = self.frames[cont]
		frame.tkraise()

# first window frame startpage
        
        
    

class NavBar(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='secondary')
        self.pack(side=LEFT, fill=Y)

        system01 = ttk.Button(self, style='info.Solid.TButton', text='Main', width=15,
                             command=lambda: Main_frame(frameprincipal))
        system01.pack(side=TOP, padx=10, pady=[250, 10])

        system02 = ttk.Button(self, style='info.Solid.TButton', text='Propeller', width=15,
                              command=lambda: Propeller_frame(frameprincipal))
        system02.pack(side=TOP, pady=[0, 10])



  


if __name__ == "__main__":

    app = ttk.Window("PAPPL", "cyborg", resizable=(False, False), size=[1080,720], iconphoto='icon.png')
    app.place_window_center()

    style = Style()
    style.configure('main.TLabel', font=('Nunito Sans', 40))
    
    MainProgram(app)

    app.mainloop()
