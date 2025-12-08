# airship_gui.py (Updated to include LOTTE Bi-Lobe shape)

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import os
import geometry_handler # Assumes geometry_handler.py is in the same directory

# --- GLOBAL CONSTANTS ---
# NOTE: The LOTTE geometry uses Gertler parameters but has the "gertler_bilobe"
# type to trigger the robust bi-lobe construction logic in geometry_handler.py.
SHAPE_DEFINITIONS = {
    "GNVR (L/D=3.044)": {"m": 0.415, "r0": 0.600, "r1": 0.180, "cp": 0.615, "l2d": 3.044, "type": "gertler"},
    "ZHIYUAN-1 (L/D=3.266)": {"m": 0.419, "r0": 0.337, "r1": 0.251, "cp": 0.651, "l2d": 3.266, "type": "gertler"},
    "Wang (L/D=3.859)": {"m": 0.404, "r0": 0.600, "r1": 0.100, "cp": 0.610, "l2d": 3.859, "type": "gertler"},
    "NPL (L/D=4.000)": {"m": 0.432, "r0": 0.589, "r1": 0.425, "cp": 0.667, "l2d": 4.000, "type": "gertler"},
    "Sphere (L/D=1.000)": {"m": 0.500, "r0": 0.500, "r1": 0.500, "cp": 0.667, "l2d": 1.000, "type": "gertler"},
    # VITAL: This line adds the LOTTE option to the GUI dropdown
    "LOTTE (L/D=3.902)": {"m": 0.4502, "r0": 0.5759, "r1": 0.1000, "cp": 0.5170, "l2d": 3.902, "type": "gertler_bilobe"},
}


# --- AIRSHIP GENERATOR APP CLASS ---
class AirshipGeneratorApp:

    COLOR_BG_DARK = '#2c3e50'
    COLOR_FRAME_DARK = '#34495e'
    COLOR_TEXT_LIGHT = '#ecf0f1'
    COLOR_ACCENT_BLUE = '#3498db'
    COLOR_BUTTON_HOVER = '#2980b9'
    COLOR_EXIT_BUTTON = '#e74c3c'
    COLOR_EXIT_HOVER = '#c0392b'
    COLOR_PARAM_KEY = '#b3cde0'

    def __init__(self, master):
        self.master = master
        master.title("SALOME Airship Geometry Generator")
        master.configure(bg=self.COLOR_BG_DARK)
        self.is_fullscreen = True

        # --- TRUE FULLSCREEN IMPLEMENTATION ---
        screen_width = self.master.winfo_screenwidth()
        screen_height = self.master.winfo_screenheight()

        self.master.geometry(f'{screen_width}x{screen_height}+0+0')
        self.master.overrideredirect(True)
        # ------------------------------------

        self.master.bind('<Escape>', lambda event: self._toggle_fullscreen())

        # --- Style Setup ---
        self.style = ttk.Style(master)
        self.style.theme_use('clam')
        DEFAULT_FONT = ('Helvetica', 14)
        TITLE_FONT = ('Roboto', 30, 'bold')
        LABEL_HEADER_FONT = ('Helvetica', 18, 'bold')
        PARAM_FONT = ('Courier', 14, 'bold')
        self.PARAM_FONT = PARAM_FONT

        # Apply styles
        self.style.configure('.', font=DEFAULT_FONT, background=self.COLOR_BG_DARK, foreground=self.COLOR_TEXT_LIGHT)
        self.style.configure('TFrame', background=self.COLOR_BG_DARK)
        self.style.configure('Card.TFrame', background=self.COLOR_FRAME_DARK, borderwidth=1, relief='solid')
        self.style.configure('TLabel', background=self.COLOR_BG_DARK, foreground=self.COLOR_TEXT_LIGHT, padding=5)
        self.style.configure('Title.TLabel', font=TITLE_FONT, foreground=self.COLOR_ACCENT_BLUE)
        self.style.configure('DropdownHeader.TLabel', font=LABEL_HEADER_FONT, foreground=self.COLOR_TEXT_LIGHT)
        self.style.configure('ParamKey.TLabel', background=self.COLOR_FRAME_DARK, foreground=self.COLOR_PARAM_KEY)
        self.style.configure('Length.TEntry', fieldbackground=self.COLOR_TEXT_LIGHT, foreground='black', font=('Courier', 14))
        self.style.configure('Primary.TButton', font=DEFAULT_FONT, foreground=self.COLOR_TEXT_LIGHT, background=self.COLOR_ACCENT_BLUE, padding=[30, 12], borderwidth=0, relief='flat')
        self.style.map('Primary.TButton', background=[('active', self.COLOR_BUTTON_HOVER)], foreground=[('active', self.COLOR_TEXT_LIGHT)])
        self.style.configure('Exit.TButton', font=DEFAULT_FONT, foreground=self.COLOR_TEXT_LIGHT, background=self.COLOR_EXIT_BUTTON, padding=[30, 12], borderwidth=0, relief='flat')
        self.style.map('Exit.TButton', background=[('active', self.COLOR_EXIT_HOVER)], foreground=[('active', self.COLOR_TEXT_LIGHT)])
        self.style.configure('Secondary.TButton', font=DEFAULT_FONT, foreground=self.COLOR_TEXT_LIGHT, background=self.COLOR_FRAME_DARK, padding=[15, 8], borderwidth=0, relief='flat')
        self.style.map('Secondary.TButton', background=[('active', self.COLOR_BUTTON_HOVER)], foreground=[('active', self.COLOR_TEXT_LIGHT)])
        self.style.configure('TCombobox', font=('Helvetica', 16), fieldbackground=self.COLOR_TEXT_LIGHT, foreground='black', padding=5)
        self.master.option_add('*TCombobox*Listbox.font', DEFAULT_FONT)
        self.master.option_add('*TCombobox*Listbox.background', self.COLOR_FRAME_DARK)
        self.master.option_add('*TCombobox*Listbox.foreground', self.COLOR_TEXT_LIGHT)

        self.selected_shape = tk.StringVar(master)
        self.selected_shape.set(list(SHAPE_DEFINITIONS.keys())[0])
        self.hull_length_var = tk.StringVar(master)
        self.hull_length_var.set("200.0")
        self.num_petals_var = tk.StringVar(master)
        self.num_petals_var.set("16")
        self.dat_file_path_var = tk.StringVar(master)
        self.dat_file_path_var.set("Click 'Browse' to set save location")


        self._create_widgets()

    def _toggle_fullscreen(self, event=None):
        self.is_fullscreen = not self.is_fullscreen

        if self.is_fullscreen:
            # Full screen mode (removes decorations, fills screen)
            screen_width = self.master.winfo_screenwidth()
            screen_height = self.master.winfo_screenheight()

            self.master.geometry(f'{screen_width}x{screen_height}+0+0')
            self.master.overrideredirect(True)
        else:
            # Windowed mode (restores decorations)
            self.master.overrideredirect(False)
            self.master.geometry("1000x800")


    def _browse_save_location(self):
        shape_name = self.selected_shape.get()
        # Create an initial filename based on current inputs
        safe_name = shape_name.replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_")
        try:
            hull_length_int = int(float(self.hull_length_var.get().strip()))
        except ValueError:
            hull_length_int = 200 # Default if input is bad

        initial_file_name = f"{safe_name}_Petal_L{hull_length_int}.dat"

        # Launch File Save As dialog
        f = filedialog.asksaveasfilename(
            defaultextension=".dat",
            initialfile=initial_file_name,
            filetypes=[("DAT files", "*.dat"), ("All files", "*.*")],
            title="Save Airship Petal Coordinates File"
        )
        if f:
            self.dat_file_path_var.set(f)


    def _create_widgets(self):
        main_container = ttk.Frame(self.master, padding="60")
        main_container.pack(expand=True, fill='both')

        main_container.grid_columnconfigure(0, weight=1)

        row_index = 0
        ttk.Label(main_container, text="🚀 Airship Geometry Generator", style='Title.TLabel').grid(row=row_index, column=0, pady=(0, 40), sticky='n')
        row_index += 1

        # --- Shape Selection ---
        ttk.Label(main_container, text="Select Standard Hull Shape:", style='DropdownHeader.TLabel').grid(row=row_index, column=0, pady=(10, 5), sticky='n')
        row_index += 1
        shape_options = list(SHAPE_DEFINITIONS.keys())
        self.shape_menu = ttk.Combobox(main_container, textvariable=self.selected_shape, values=shape_options, state="readonly", width=30)
        self.shape_menu.grid(row=row_index, column=0, pady=(0, 20), padx=10, sticky='n')
        self.shape_menu.bind("<<ComboboxSelected>>", self._show_parameters)
        row_index += 1

        # --- Hull and Petal Input Frame ---
        input_frame = ttk.Frame(main_container, style='Card.TFrame', padding="18 12 18 12")
        input_frame.grid(row=row_index, column=0, pady=(10, 12), padx=18, sticky='ew')
        input_frame.grid_columnconfigure(0, weight=1)
        input_frame.grid_columnconfigure(1, weight=1)

        ttk.Label(input_frame, text="Hull Length [units]:", font=('Helvetica', 14, 'bold'), background=self.COLOR_FRAME_DARK, foreground=self.COLOR_ACCENT_BLUE).grid(row=0, column=0, sticky='w', pady=5, padx=10)
        self.hull_length_entry = ttk.Entry(input_frame, textvariable=self.hull_length_var, style='Length.TEntry', width=10)
        self.hull_length_entry.grid(row=0, column=1, sticky='e', pady=5, padx=10)

        ttk.Label(input_frame, text="Number of Hull Petals:", font=('Helvetica', 14, 'bold'), background=self.COLOR_FRAME_DARK, foreground=self.COLOR_ACCENT_BLUE).grid(row=1, column=0, sticky='w', pady=5, padx=10)
        self.num_petals_entry = ttk.Entry(input_frame, textvariable=self.num_petals_var, style='Length.TEntry', width=10)
        self.num_petals_entry.grid(row=1, column=1, sticky='e', pady=5, padx=10)

        row_index += 1

        # --- File Save Frame ---
        file_frame = ttk.Frame(main_container, style='Card.TFrame', padding="18 12 18 12")
        file_frame.grid(row=row_index, column=0, pady=(10, 12), padx=18, sticky='ew')
        file_frame.grid_columnconfigure(0, weight=1)
        file_frame.grid_columnconfigure(1, weight=3)
        file_frame.grid_columnconfigure(2, weight=1)

        ttk.Label(file_frame, text="DAT Save Location:", font=('Helvetica', 14, 'bold'), background=self.COLOR_FRAME_DARK, foreground=self.COLOR_ACCENT_BLUE).grid(row=0, column=0, sticky='w', padx=10)

        ttk.Label(file_frame, textvariable=self.dat_file_path_var, font=('Courier', 10), background=self.COLOR_FRAME_DARK, foreground=self.COLOR_TEXT_LIGHT, wraplength=300).grid(row=0, column=1, sticky='ew', padx=5)

        ttk.Button(file_frame, text="Browse", command=self._browse_save_location, style='Secondary.TButton').grid(row=0, column=2, sticky='e', padx=10)
        row_index += 1

        # --- Design Variables Frame ---
        self.param_frame = ttk.Frame(main_container, padding="30 25 30 25", style='Card.TFrame')
        self.param_frame.grid(row=row_index, column=0, pady=40, padx=20, sticky='ew')
        self.param_labels = {}
        row_index += 1

        self._show_parameters(None)

        # --- Action Buttons Frame ---
        button_frame = ttk.Frame(main_container)
        button_frame.grid(row=row_index, column=0, pady=(20, 40), padx=20, sticky='ew')
        button_frame.grid_columnconfigure(0, weight=1)
        button_frame.grid_columnconfigure(1, weight=1)

        ttk.Button(button_frame, text="Generate Geometry & DAT File", command=self._generate_and_run, style='Primary.TButton').grid(row=0, column=0, padx=25, sticky='ew')
        ttk.Button(button_frame, text="Exit Application", command=self.master.destroy, style='Exit.TButton').grid(row=0, column=1, padx=25, sticky='ew')
        row_index += 1

        # --- Fullscreen Toggle Button (outside main_container) ---
        ttk.Button(self.master, text="Toggle Fullscreen (ESC)", command=self._toggle_fullscreen, style='Secondary.TButton').pack(side='bottom', pady=(0, 20))


    def _show_parameters(self, event):
        # Clears and redraws the Design Variables based on the selected shape
        for widget in self.param_frame.winfo_children():
            widget.destroy()

        shape_name = self.selected_shape.get()
        params = SHAPE_DEFINITIONS.get(shape_name, {})

        ttk.Label(self.param_frame, text="Design Variables:", font=('Helvetica', 18, 'bold'), background=self.COLOR_FRAME_DARK, foreground=self.COLOR_ACCENT_BLUE).grid(row=0, column=0, columnspan=2, pady=(0, 20), sticky='w')

        row = 1
        # Filter out the 'type' key before display, as it's only for backend logic
        display_params = {k: v for k, v in params.items() if k != 'type'}

        for key, value in display_params.items():
            ttk.Label(self.param_frame, text=f"{key.upper()}:", style='ParamKey.TLabel', font=('Helvetica', 14)).grid(row=row, column=0, sticky='e', padx=(0, 15), pady=3)
            ttk.Label(self.param_frame, text=f"{value:.4f}" if isinstance(value, float) else f"{value}", font=self.PARAM_FONT, background=self.COLOR_FRAME_DARK).grid(row=row, column=1, sticky='w', pady=3)
            row += 1

        self.param_frame.grid_columnconfigure(0, weight=1, minsize=150)
        self.param_frame.grid_columnconfigure(1, weight=1, minsize=100)


    def _validate_inputs(self):
        # 1. Validate Hull Length
        hull_length_str = self.hull_length_var.get().strip()
        try:
            hull_length = float(hull_length_str)
            if hull_length <= 0:
                raise ValueError("Hull length must be a positive number.")
        except Exception:
            messagebox.showerror("Input Error", "Please enter a valid positive number for Hull Length.")
            return False, 0, 0, ""

        # 2. Validate Number of Petals
        num_petals_str = self.num_petals_var.get().strip()
        try:
            num_petals = int(num_petals_str)
            if num_petals < 4 or num_petals > 100: # Setting a reasonable upper limit
                raise ValueError("Number of petals must be an integer between 4 and 100.")
        except Exception:
            messagebox.showerror("Input Error", "Please enter a valid integer (4-100) for Number of Petals.")
            return False, 0, 0, ""

        # 3. Validate File Path
        dat_filename = self.dat_file_path_var.get().strip()
        if not dat_filename or dat_filename == "Click 'Browse' to set save location":
            messagebox.showerror("Input Error", "Please use the 'Browse' button to select a valid file save location.")
            return False, 0, 0, ""

        return True, hull_length, num_petals, dat_filename

    # --- Execution function calling the handler ---
    def _generate_and_run(self):

        # 1. Validate all inputs
        is_valid, hull_length, num_petals, dat_filename = self._validate_inputs()
        if not is_valid:
            return

        shape_name = self.selected_shape.get()
        # Retrieve the parameters dictionary
        selected_params = SHAPE_DEFINITIONS.get(shape_name)

        # 2. Call the external function from geometry_handler.py
        result_message = geometry_handler.generate_and_run_geometry(
            shape_name,
            selected_params,
            hull_length,
            num_petals,
            dat_filename
        )

        # 3. Display the combined result from the handler
        messagebox.showinfo("Processing Results", result_message)


if __name__ == "__main__":
    try:
        root = tk.Tk()
        # Ensure correct DPI scaling on Windows
        try:
            from ctypes import windll
            windll.shcore.SetProcessDpiAwareness(1)
        except:
            pass

        app = AirshipGeneratorApp(root)
        root.mainloop()
    except ImportError as e:
        messagebox.showerror("Import Error", f"A required library (tkinter, numpy, or subprocess) is missing. Error: {e}")