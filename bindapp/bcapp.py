import tkinter
import tkinter.messagebox
import customtkinter
import webbrowser
from tkinter import filedialog
import subprocess
import pickle
import os
from PIL import Image, ImageTk

# Modes: "System" (standard), "Dark", "Light"
customtkinter.set_appearance_mode("System")

# Themes: "blue" (standard), "green", "dark-blue"
customtkinter.set_default_color_theme("blue")


class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        with open("bindapp/help.p", "rb") as handle:
            self.help = pickle.load(handle)

        # Configure Heading
        self.title("BindCompare")
        self.geometry(f"{1100}x{600}")

        # Configure the Layout
        # configure grid layout (4x4)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)
        self.grid_rowconfigure((0, 1), weight=1)

        # create sidebar frame with widgets
        self.sidebar_frame = customtkinter.CTkFrame(self, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
        self.sidebar_frame.grid_rowconfigure(4, weight=1)
        self.logo_label = customtkinter.CTkLabel(
            self.sidebar_frame,
            text="BindCompare",
            font=customtkinter.CTkFont(size=20, weight="bold"),
        )
        self.logo_label.grid(row=0, column=0, padx=20, pady=(20, 10))
        self.sidebar_button_1 = customtkinter.CTkButton(
            self.sidebar_frame, command=self.launch_manual, text="Manual Page"
        )
        self.sidebar_button_1.grid(row=3, column=0, padx=20, pady=10)
        self.appearance_mode_label = customtkinter.CTkLabel(
            self.sidebar_frame, text="Appearance Mode:", anchor="w"
        )
        self.appearance_mode_label.grid(row=5, column=0, padx=20, pady=(10, 0))
        self.appearance_mode_optionemenu = customtkinter.CTkOptionMenu(
            self.sidebar_frame,
            values=["Light", "Dark", "System"],
            command=self.change_appearance_mode_event,
        )
        self.appearance_mode_optionemenu.grid(row=6, column=0, padx=20, pady=(10, 10))
        self.scaling_label = customtkinter.CTkLabel(
            self.sidebar_frame, text="UI Scaling:", anchor="w"
        )
        self.scaling_label.grid(row=7, column=0, padx=20, pady=(10, 0))
        self.scaling_optionemenu = customtkinter.CTkOptionMenu(
            self.sidebar_frame,
            values=["80%", "90%", "100%", "110%", "120%"],
            command=self.change_scaling_event,
        )
        self.scaling_optionemenu.grid(row=8, column=0, padx=20, pady=(10, 20))

        # MAIN PAGE
        ## INPUTS
        self.input_frame = customtkinter.CTkScrollableFrame(
            self, label_text="Enter BindCompare Inputs"
        )
        self.input_frame.grid(
            row=0, column=1, padx=(20, 20), pady=(20, 20), sticky="nsew"
        )
        #### add all input options
        self.ref_entry = customtkinter.CTkEntry(
            self.input_frame, placeholder_text="Reference BED Filepath", width=300
        )
        self.ref_entry.grid(row=2, column=0, padx=10, pady=5, sticky="nsew")
        self.ref_button = customtkinter.CTkButton(
            self.input_frame,
            text="Choose File",
            command=lambda entry=self.ref_entry: self.open_file(entry),
        )
        self.ref_button.grid(row=2, column=1, padx=10, pady=5)
        #### exp entry
        self.exp_entry = customtkinter.CTkEntry(
            self.input_frame, placeholder_text="Overlayed BED Filepath", width=300
        )
        self.exp_entry.grid(row=4, column=0, padx=10, pady=5, sticky="nsew")
        self.exp_button = customtkinter.CTkButton(
            self.input_frame,
            text="Choose File",
            command=lambda entry=self.exp_entry: self.open_file(entry),
        )
        self.exp_button.grid(row=4, column=1, padx=10, pady=5)
        #### scope entry
        vcmd = self.input_frame.register(self.intcallback)
        self.scope = customtkinter.CTkEntry(
            self.input_frame, placeholder_text="Enter Numerical Scope", width=300
        )
        self.scope.configure(validate="key", validatecommand=(vcmd, "%P"))
        self.scope_button = customtkinter.CTkButton(
            self.input_frame, text="Choose Scope", state="disabled"
        )
        self.scope.grid(row=6, column=0, padx=10, pady=5, sticky="nsew")
        self.scope_button.grid(row=6, column=1, padx=10, pady=5, sticky="nsew")
        #### name entry
        self.name = customtkinter.CTkEntry(
            self.input_frame, placeholder_text="BCExpName", width=300
        )
        self.name_button = customtkinter.CTkButton(
            self.input_frame, text="Choose Name", state="disabled"
        )
        self.name.grid(row=8, column=0, padx=10, pady=5, sticky="nsew")
        self.name_button.grid(row=8, column=1, padx=10, pady=5, sticky="nsew")
        #### out-dir
        self.outdir_entry = customtkinter.CTkEntry(
            self.input_frame, placeholder_text="Output Directory", width=300
        )
        self.outdir_entry.grid(row=10, column=0, padx=10, pady=5, sticky="nsew")
        self.out_button = customtkinter.CTkButton(
            self.input_frame,
            text="Choose Folder",
            command=lambda entry=self.outdir_entry: self.open_file(entry, True),
        )
        self.out_button.grid(row=10, column=1, padx=10, pady=5)
        #### Gene GTF File
        self.gtf_entry = customtkinter.CTkEntry(
            self.input_frame, placeholder_text="[Optional] Genes GTF File", width=300
        )
        self.gtf_entry.grid(row=12, column=0, padx=10, pady=5, sticky="nsew")
        self.gtf_button = customtkinter.CTkButton(
            self.input_frame,
            text="Choose File",
            command=lambda entry=self.gtf_entry: self.open_file(entry),
        )
        self.gtf_button.grid(row=12, column=1, padx=10, pady=5)
        #### FASTA File
        self.fa_entry = customtkinter.CTkEntry(
            self.input_frame, placeholder_text="[Optional] Genome FA File", width=300
        )
        self.fa_entry.grid(row=14, column=0, padx=10, pady=5, sticky="nsew")
        self.fa_button = customtkinter.CTkButton(
            self.input_frame,
            text="Choose File",
            command=lambda entry=self.fa_entry: self.open_file(entry),
        )
        self.fa_button.grid(row=14, column=1, padx=10, pady=5, sticky="nsew")

        ## Add Run BindCompare
        self.run_frame = customtkinter.CTkFrame(self)
        self.run_frame.grid(row=1, column=1, padx=20, pady=(0, 20), sticky="nsew")
        self.run_frame.grid_columnconfigure(2, weight=1)
        self.run_frame.grid_rowconfigure((0, 1), weight=1)

        self.run_entry = customtkinter.CTkEntry(
            self.run_frame, placeholder_text="BindCompare Currently IDLE.", width=300
        )
        self.run_button = customtkinter.CTkButton(
            self.run_frame,
            text="Run BindCompare",
            command=lambda entry=self.run_entry: self.launch_analysis(entry),
        )
        self.run_entry.grid(row=0, column=1, padx=20, pady=0)
        self.run_button.grid(row=0, column=0, padx=20, pady=0)
        ## Add Text Entry with Manual
        self.help_manual = customtkinter.CTkTextbox(self.run_frame, height=100)
        self.help_manual.grid(
            row=1, column=0, columnspan=2, padx=20, pady=(0, 20), sticky="nsew"
        )
        self.help_manual.insert(tkinter.END, self.help["General Help"])
        self.help_manual.configure(state="disabled")

        ## Add Visualizer Panel
        self.visualize_frame = customtkinter.CTkScrollableFrame(
            self, label_text="Visualize and Interpret Results", width=300
        )
        self.visualize_frame.grid(
            row=0, rowspan=2, column=2, padx=(0, 20), pady=20, sticky="nsew"
        )
        self.visualize_frame.grid_columnconfigure(0, weight=0)
        self.visualize_frame.grid_columnconfigure((1, 2), weight=1)
        self.directory_label = customtkinter.CTkLabel(
            self.visualize_frame, text="Out. Dir"
        )
        self.directory_label.grid(
            row=0, column=0, padx=(10, 10), pady=20, sticky="nsew"
        )
        self.directory_entry = customtkinter.CTkEntry(self.visualize_frame)
        self.directory_entry.grid(row=0, column=1, padx=0, pady=20, sticky="nsew")
        self.browse_dir = customtkinter.CTkButton(
            self.visualize_frame, text="Browse", command=self.choose_directory
        )
        self.browse_dir.grid(row=0, column=2, padx=(10, 10), pady=20)
        self.image_dropdown_label = customtkinter.CTkLabel(
            self.visualize_frame, text="Select Image:"
        )
        self.image_dropdown_label.grid(row=1, column=0, columnspan=2, padx=(0, 50))

        self.image_dropdown = customtkinter.CTkComboBox(self.visualize_frame, width=200)
        self.image_dropdown.grid(row=1, column=1, columnspan=2, padx=(44, 0))

        self.show_image_button = customtkinter.CTkButton(
            self.visualize_frame, text="Show Image", command=self.show_selected_image
        )
        self.show_image_button.grid(
            row=2, column=0, padx=(20, 10), pady=(10, 20), columnspan=3
        )

        self.image_label = customtkinter.CTkLabel(self.visualize_frame, text="")
        self.image_label.grid(row=3, column=0, columnspan=3, padx=20, pady=10)

        self.interpret_label = customtkinter.CTkLabel(
            self.visualize_frame, text="Select Plot Interpretation Category"
        )
        self.interpret_label.grid(row=4, column=0, columnspan=3)

        help_options = [key for key in self.help.keys()]
        help_options.remove("General Help")
        self.interpret_menu = customtkinter.CTkComboBox(
            self.visualize_frame,
        )
        self.interpret_menu.configure(values=help_options)
        self.interpret_menu.grid(
            row=5, column=0, columnspan=2, padx=(0, 30), pady=(20, 0)
        )
        self.interpret_run = customtkinter.CTkButton(
            self.visualize_frame, text="Show Help Text", command=self.show_help
        )
        self.interpret_run.grid(
            row=5, column=1, columnspan=2, padx=(60, 0), pady=(20, 0)
        )

        self.int_help_box = customtkinter.CTkTextbox(self.visualize_frame, width=250)
        self.int_help_box.grid(
            row=6, column=0, columnspan=3, pady=20, padx=20, sticky="nsew"
        )

        ## Jaccard Similarity

    def launch_manual(self):
        """
        Launch BindCompare README for manual page!
        """
        webbrowser.open("https://github.com/pranavmahabs/bindcompare")

    def open_file(self, entry, dir=False):
        """
        Read in the file path
        """
        if dir:
            file_path = filedialog.askdirectory()
        else:
            file_path = filedialog.askopenfilename()
        if file_path:
            # filename = file_path.split("/")[-1]  # Extracting just the filename
            entry.delete(0, tkinter.END)
            entry.insert(tkinter.END, file_path)

    def change_appearance_mode_event(self, new_appearance_mode: str):
        customtkinter.set_appearance_mode(new_appearance_mode)

    def change_scaling_event(self, new_scaling: str):
        new_scaling_float = int(new_scaling.replace("%", "")) / 100
        customtkinter.set_widget_scaling(new_scaling_float)

    def intcallback(self, P):
        if str.isdigit(P) or P == "":
            return True
        else:
            return False

    def launch_analysis(self, entry):
        entry.delete(0, tkinter.END)
        entry.insert(tkinter.END, "Launching BindCompare... status in terminal.")

        # Launching Merge
        subprocess.run(
            [
                "./bindapp/bindcompare.sh",
                self.ref_entry.get(),
                self.exp_entry.get(),
                self.scope.get(),
                self.name.get(),
                self.outdir_entry.get(),
                self.gtf_entry.get(),
                self.fa_entry.get(),
            ]
        )
        entry.delete(0, tkinter.END)
        entry.insert(tkinter.END, "BindCompare COMPLETED.")

    def choose_directory(self):
        directory_path = filedialog.askdirectory()
        self.directory_entry.delete(0, tkinter.END)
        self.directory_entry.insert(0, directory_path)
        self.populate_dropdown(directory_path)

    def populate_dropdown(self, directory_path):
        image_files = [
            file
            for file in os.listdir(directory_path)
            if file.endswith(("png", "jpg", "jpeg", "gif"))
        ]
        # self.image_dropdown["values"] = image_files
        self.image_dropdown.configure(values=image_files)

    def show_selected_image(self):
        selected_image = self.image_dropdown.get()
        if selected_image:
            image_path = os.path.join(self.directory_entry.get(), selected_image)
            self.load_and_display_image(image_path)

    def load_and_display_image(self, image_path):
        image = customtkinter.CTkImage(
            light_image=Image.open(image_path), size=(250, 250)
        )
        self.image_label.configure(image=image)
        # self.image_label.image = photo  # Keep a reference to prevent garbage collection

    def show_help(self):
        selected_option = self.interpret_menu.get()
        help_text = self.help.get(selected_option, "No help text available.")
        self.int_help_box.delete(1.0, tkinter.END)
        self.int_help_box.insert(tkinter.END, help_text)


if __name__ == "__main__":
    app = App()
    app.mainloop()