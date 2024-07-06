import tkinter
import tkinter.messagebox
import customtkinter
import webbrowser
from tkinter import filedialog
import subprocess
import pickle
import os
from PIL import Image, ImageTk

from .bindapp.help_dict import HELP_DICT
from .comparexp import verify_summary_file

# Modes: "System" (standard), "Dark", "Light"
customtkinter.set_appearance_mode("System")

# Themes: "blue" (standard), "green", "dark-blue"
customtkinter.set_default_color_theme("blue")


class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        # with open("bindapp/help.p", "rb") as handle:
        #     self.help = pickle.load(handle)
        self.help = HELP_DICT

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

        self.sidebar_button_2 = customtkinter.CTkButton(
            self.sidebar_frame,
            command=self.open_new_window,
            text="comparexp",
        )
        self.sidebar_button_2.grid(row=4, column=0, padx=20, pady=10)

        self.appearance_mode_label = customtkinter.CTkLabel(
            self.sidebar_frame, text="Appearance Mode:", anchor="w"
        )
        self.appearance_mode_label.grid(row=5, column=0, padx=20, pady=(10, 0))
        self.appearance_mode_optionemenu = customtkinter.CTkOptionMenu(
            self.sidebar_frame,
            values=["System", "Dark", "Light"],
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

        self.image_dropdown = customtkinter.CTkComboBox(
            self.visualize_frame,
            width=200,
            values=["Image Dropdown"],
        )
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
            values=["Help Menu"],
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

    def open_new_window(self):
        new_window = customtkinter.CTkToplevel(self)  # Creating a new window
        new_window.title("comparexp")  # Setting title for the new window

        new_window.grid_columnconfigure(0, weight=1)
        new_window.grid_rowconfigure(0, weight=1)

        analysis_frame = customtkinter.CTkFrame(new_window, corner_radius=0)
        analysis_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        analysis_frame.grid_columnconfigure((0, 1), weight=1)
        analysis_frame.grid_rowconfigure((5), weight=1)

        self.entry_one = customtkinter.CTkEntry(
            analysis_frame, placeholder_text="BC Output Directory Path 1", width=300
        )
        self.entry_one.grid(row=0, column=0, padx=10, pady=5, sticky="nsew")

        self.entry_one_button = customtkinter.CTkButton(
            analysis_frame,
            text="BC Output Dir 1",
            command=self.choose_bcdirectory1,
        )
        self.entry_one_button.grid(row=0, column=1, padx=10, pady=5, sticky="nsew")

        self.entry_two = customtkinter.CTkEntry(
            analysis_frame, placeholder_text="BC Output Directory Path 2", width=300
        )
        self.entry_two.grid(row=1, column=0, padx=10, pady=5, sticky="nsew")

        self.entry_two_button = customtkinter.CTkButton(
            analysis_frame,
            text="BC Output Dir 2",
            command=self.choose_bcdirectory2,
        )
        self.entry_two_button.grid(row=1, column=1, padx=10, pady=5, sticky="nsew")

        self.comparexp = customtkinter.CTkButton(
            analysis_frame,
            text="Run comparexp",
            command=self.launch_comparexp,
        )
        self.comparexp.grid(
            row=2, column=0, padx=10, pady=10, columnspan=2, sticky="nsew"
        )

        self.image_label2 = customtkinter.CTkLabel(analysis_frame, text="")
        self.image_label2.grid(row=3, column=0, columnspan=2, padx=20, pady=10)

        self.label_comparexp = customtkinter.CTkLabel(
            analysis_frame,
            text="Errors Will be Shown in Terminal/Command Line. Open only one comparexp at a time.",
        )
        self.label_comparexp.grid(row=4, column=0, columnspan=2, padx=10, pady=10)

        self.compare_out = customtkinter.CTkTextbox(analysis_frame, height=40)
        self.compare_out.grid(
            row=5, column=0, columnspan=2, padx=20, pady=10, sticky="nsew"
        )

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

        gtf = self.gtf_entry.get()
        if len(self.gtf_entry.get()) == 0:
            gtf = "None"
        fa = self.fa_entry.get()
        if len(self.fa_entry.get()) == 0:
            fa = "None"

        # Launching Merge
        subprocess.run(
            [
                "bindcompare",
                "--ref",
                self.ref_entry.get(),
                "--exp",
                self.exp_entry.get(),
                "--scope",
                self.scope.get(),
                "--name",
                self.name.get(),
                "--out",
                self.outdir_entry.get(),
                "--gtf",
                gtf,
                "--fasta",
                fa,
            ]
        )
        entry.delete(0, tkinter.END)
        entry.insert(tkinter.END, "BindCompare COMPLETED.")

    def launch_comparexp(self):
        subprocess.run(
            [
                "comparexp",
                "--bindpath_1",
                self.entry_one.get(),
                "--bindpath_2",
                self.entry_two.get(),
            ]
        )
        self.show_selected_image2()
        self.fill_comparexp()

    def choose_directory(self):
        directory_path = filedialog.askdirectory()
        self.directory_entry.delete(0, tkinter.END)
        self.directory_entry.insert(0, directory_path)
        self.populate_dropdown(directory_path)

    def choose_bcdirectory1(self):
        directory_path = filedialog.askdirectory()
        self.entry_one.delete(0, tkinter.END)
        self.entry_one.insert(0, directory_path)

    def choose_bcdirectory2(self):
        directory_path = filedialog.askdirectory()
        self.entry_two.delete(0, tkinter.END)
        self.entry_two.insert(0, directory_path)

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

    def get_comparexp_name(self):
        summary_files_folder1 = verify_summary_file(self.entry_one.get())
        summary_files_folder2 = verify_summary_file(self.entry_two.get())

        if summary_files_folder1 and summary_files_folder2:
            p1 = os.path.basename(summary_files_folder1[0]).split("_summary.txt")[0]
            p2 = os.path.basename(summary_files_folder2[0]).split("_summary.txt")[0]
            return f"{p1}_v_{p2}_venn.png", f"{p1}_v_{p2}_summary.txt"

    def show_selected_image2(self):
        image_path = self.get_comparexp_name()[0]
        if not os.path.exists(image_path):
            print(f"'{image_path}' not found. comparexp likely failed.")
        else:
            self.load_and_display_image2(image_path)

    def load_and_display_image2(self, image_path):
        image = customtkinter.CTkImage(
            light_image=Image.open(image_path), size=(420, 280)
        )
        self.image_label2.configure(image=image)
        # self.image_label.image = photo  # Keep a reference to prevent garbage collection

    def fill_comparexp(self):
        text_path = self.get_comparexp_name()[1]
        text = ""
        if not os.path.exists(text_path):
            text = f"'{text_path}' not found. comparexp likely failed."
        else:
            with open(text_path, "r") as handle:
                for line in handle:
                    text += line + "\n"
        self.compare_out.delete(1.0, tkinter.END)
        self.compare_out.insert(tkinter.END, text)

    def show_help(self):
        selected_option = self.interpret_menu.get()
        help_text = self.help.get(selected_option, "No help text available.")
        self.int_help_box.delete(1.0, tkinter.END)
        self.int_help_box.insert(tkinter.END, help_text)


def main():
    app = App()
    app.mainloop()


if __name__ == "__main__":
    app = App()
    app.mainloop()
