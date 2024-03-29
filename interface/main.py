# Import necessary libraries
import tkinter as tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from ttkbootstrap import Style
from main_page import main_page
from main_wthout import main_wthout

# Define the main program class
class MainProgram(ttk.Frame):
    def __init__(self, master_window):
        super().__init__(master_window)
        self.pack(fill=BOTH, expand=YES)
        
        # Create navigation bar and dowbar (database bar)
        NavBar(self)
        self.db = DowBar(self)
        
        # Create a list of pages depending on the number of systems in the navigation bar
        self.pages = [main_page(self), main_wthout(self)]
        self.actual_page = 0
        
        # Show dowbar only if the actual page is not the main page
        if(self.actual_page != 0):
            self.db.show()

        # Show the initial page
        self.show_pages(0)
        self.values = []

    # Switch between pages
    def show_pages(self, page):
        self.pages[self.actual_page].unshow()
        self.actual_page = page
        frame = self.pages[self.actual_page]
        
        # Show dowbar only if the page is not the main page
        if(page != 0):
            self.db.show()
        else:
            self.db.unshow()
            
        frame.show()

    # Show values from the main page
    def executer(self):
        self.pages[0].executer_program()

# Define the DowBar class
class DowBar(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='secondary')
        
        frame = ttk.Frame(self, style='secondary')
        frame.pack()

    # Show the DowBar
    def show(self):
        self.pack(side=BOTTOM, fill=BOTH)
        
    # Hide the DowBar
    def unshow(self):
        self.pack_forget()

# Define the Navigation Bar class
class NavBar(ttk.Frame):
    def __init__(self, frameprincipal):
        super().__init__(frameprincipal, style='secondary')
        self.pack(side=LEFT,fill=Y)
        
        # Create buttons for navigating to different pages
        system01 = ttk.Button(self, style='info.Solid.TButton', text='Main',width=15, command = lambda : frameprincipal.show_pages(0))
        system01.pack(side=TOP,padx=10,pady=[250,10])

        system02 = ttk.Button(self, style='info.Solid.TButton', text='Main wthout',width=15, command = lambda : frameprincipal.show_pages(1))
        system02.pack(side=TOP,padx=10,pady=[0,10])

        

# Create the main application window
app = ttk.Window("PAPPL", "cyborg", resizable=(False, False), size=[1080,720], iconphoto='icon.png')
app.place_window_center()

# Create a style for the application
style = Style()
style.configure('main.TLabel', font=('Arial', 40))

# Initialize the main program
main_program_instance = MainProgram(app)

# Start the application event loop
app.mainloop()
