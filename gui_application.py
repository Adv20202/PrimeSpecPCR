import os
import sys
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import subprocess
import threading
import queue
import time
import re
import configparser
import signal
from pathlib import Path
import datetime

class PrimeSpecPCRApp:
    def restart_application(self):
        """Restart the entire application."""
        if self.running:
            if not messagebox.askyesno("Warning", "A process is currently running. Stop it and restart the application?"):
                return
            self.stop_process()
        
        self.update_log("\nRestarting application...\n")
        
        # Clean up any input files before restart
        self.cleanup_input_files()
        
        # Python executable and this script path
        python = sys.executable
        script = os.path.abspath(sys.argv[0])
        
        # Close the current application
        self.root.destroy()
        
        # Start a new instance
        os.execl(python, python, script)


def main():
    root = tk.Tk()
    
    # Set application icon if available
    try:
        icon_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "icon.ico")
        if os.path.exists(icon_path):
            root.iconbitmap(icon_path)
    except:
        pass
        
    # Configure styles for input message with yellow background
    style = ttk.Style()
    style.configure("YellowBg.TLabel", background="#FFFF99", foreground="black")
    
    app = PrimeSpecPCRApp(root)
    root.mainloop()

class PrimeSpecPCRApp:
    def __init__(self, root):
        self.root = root
        self.root.title("PrimeSpecPCR: Species-Specific Primer Design Toolkit")
        
        screen_width = root.winfo_screenwidth()
        
        self.root.geometry(f"{screen_width}x900")  
        self.root.minsize(1200, 750)
        
        # Program parameters
        self.email = tk.StringVar()
        self.api_key = tk.StringVar()
        self.taxids = tk.StringVar()
        self.load_config()
        
        # Variable to track current module
        self.current_module = 1
        self.process = None
        self.output_queue = queue.Queue()
        self.input_queue = queue.Queue()
        self.timeout_id = None
        self.last_output_time = None
        
        # Application state
        self.running = False
        self.waiting_for_input = False
        
        # For file-based input handling
        self.input_checker_id = None
        
        # Create main interface
        self.create_gui()
        
        # Show initial help message
        self.show_initial_help()
        
    def load_config(self):
        """Load saved configuration."""
        config = configparser.ConfigParser()
        config_path = Path('user_config.ini')
        
        if config_path.exists():
            config.read(config_path)
            if 'User' in config:
                self.email.set(config.get('User', 'email', fallback=''))
                self.api_key.set(config.get('User', 'api_key', fallback=''))
                self.taxids.set(config.get('User', 'taxids', fallback=''))
    
    def save_config(self):
        """Save user configuration."""
        config = configparser.ConfigParser()
        config['User'] = {
            'email': self.email.get(),
            'api_key': self.api_key.get(),
            'taxids': self.taxids.get()
        }
        
        try:
            with open('user_config.ini', 'w') as f:
                config.write(f)
            self.update_log("Configuration saved successfully.\n")
            self.show_module_1_help()
        except Exception as e:
            self.update_log(f"Error saving configuration: {e}\n")
    
    def create_gui(self):
        """Create main application interface."""
        # Create main frame
        main_frame = ttk.Frame(self.root, padding=10)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Left panel - configuration and controls
        left_frame = ttk.LabelFrame(main_frame, text="Configuration", padding=10)
        left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        
        # User parameters
        user_frame = ttk.Frame(left_frame)
        user_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(user_frame, text="Email:").pack(anchor=tk.W)
        ttk.Entry(user_frame, textvariable=self.email, width=30).pack(fill=tk.X, pady=2)
        
        ttk.Label(user_frame, text="NCBI API Key:").pack(anchor=tk.W)
        ttk.Entry(user_frame, textvariable=self.api_key, width=30).pack(fill=tk.X, pady=2)
        
        ttk.Label(user_frame, text="TaxID (comma separated):").pack(anchor=tk.W)
        ttk.Entry(user_frame, textvariable=self.taxids, width=30).pack(fill=tk.X, pady=2)
        
        ttk.Button(user_frame, text="Save Configuration", command=self.save_config).pack(anchor=tk.W, pady=5)
        
        # Separator
        ttk.Separator(left_frame, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        # Module buttons
        modules_frame = ttk.LabelFrame(left_frame, text="Program Modules", padding=10)
        modules_frame.pack(fill=tk.X, pady=5)
        
        self.module_buttons = []
        module_commands = [
            ("1. Genetic Sequence Retrieval", lambda: self.run_module(1)),
            ("2. Multiple Sequence Alignment", lambda: self.run_module(2)),
            ("3. PCR Primer Design", lambda: self.run_module(3)),
            ("4. Primer Specificity Testing", lambda: self.run_module(4))
        ]
        
        for text, command in module_commands:
            btn = ttk.Button(modules_frame, text=text, command=command)
            btn.pack(fill=tk.X, pady=2)
            self.module_buttons.append(btn)
        
        # Control buttons
        control_frame = ttk.Frame(left_frame)
        control_frame.pack(fill=tk.X, pady=10)
        
        # Stop process button
        self.stop_button = ttk.Button(
            control_frame, 
            text="Stop Current Process", 
            command=self.stop_process, 
            style="Stop.TButton"
        )
        self.stop_button.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=2)
        
        # Restart button
        self.restart_button = ttk.Button(
            control_frame,
            text="Restart Application",
            command=self.restart_application,
            style="Restart.TButton"
        )
        self.restart_button.pack(side=tk.RIGHT, fill=tk.X, expand=True, padx=2)
        
        # Add help panel to the left side
        help_frame = ttk.LabelFrame(left_frame, text="Help & Guidance", padding=10)
        help_frame.pack(fill=tk.BOTH, pady=5, expand=True)
        
        self.help_text = scrolledtext.ScrolledText(help_frame, wrap=tk.WORD, height=12, bg="#FFFFFF", fg="#000000")
        self.help_text.pack(fill=tk.BOTH, expand=True)
        
        # Right panel - logs and output
        right_frame = ttk.LabelFrame(main_frame, text="Logs and Program Output", padding=10)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Text field for logs
        self.log_text = scrolledtext.ScrolledText(right_frame, wrap=tk.WORD, height=25, bg="#FFFFFF", fg="#000000")
        self.log_text.pack(fill=tk.BOTH, expand=True)
        self.log_text.config(state=tk.DISABLED)
        
        # Create user input frame at bottom of right panel
        self.input_frame = ttk.Frame(right_frame)
        self.input_frame.pack(fill=tk.X, pady=5)
        
        self.input_label = ttk.Label(self.input_frame, text="Console input:")
        self.input_label.pack(side=tk.LEFT, padx=5)
        
        self.input_var = tk.StringVar()
        self.input_entry = ttk.Entry(self.input_frame, textvariable=self.input_var, width=30)
        self.input_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        self.input_entry.bind("<Return>", lambda event: self.submit_input())
        
        self.submit_button = ttk.Button(self.input_frame, text="Send", command=self.submit_input)
        self.submit_button.pack(side=tk.RIGHT, padx=5)
        
        # Add an input message label with yellow background for highlighting
        self.input_message = ttk.Label(right_frame, text="", foreground="black", background="#FFFF99")
        self.input_message.pack(fill=tk.X, pady=2)
        # Initially hide the input message label
        self.input_message.pack_forget()
        
        # Initially disable input controls
        self.input_entry.config(state=tk.DISABLED)
        self.submit_button.config(state=tk.DISABLED)
        
        # Status bar with progress
        status_frame = ttk.Frame(self.root)
        status_frame.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.status_var = tk.StringVar()
        self.status_var.set("Ready")
        status_bar = ttk.Label(status_frame, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        self.progress_bar = ttk.Progressbar(status_frame, mode='indeterminate', length=100)
        self.progress_bar.pack(side=tk.RIGHT, padx=5)
        
        # Style for buttons
        style = ttk.Style()
        style.configure("Stop.TButton", foreground="red")
        style.configure("Restart.TButton", foreground="blue")
        
        # Show welcome message
        self.update_log("""Welcome to PrimeSpecPCR: Species-Specific Primer Design Toolkit!

        This program consists of 4 modules that should be run sequentially. 
        Each module should be launched by clicking the button under the subtitle "Program Modules". 
        Run the modules one by one from top to bottom, but launch the next one only after the previous 
        one has finished its work — this will be indicated in this terminal window.

        To get started, enter your email, API key, and TaxID numbers, 
        then save the configuration and run the first module.

        Also, pay attention to the window under the subtitle "Help & Guidance". 
        This is a contextual help window where the content changes dynamically.
        """)

    def update_log(self, text):
        """Update log field."""
        self.log_text.config(state=tk.NORMAL)
        self.log_text.insert(tk.END, text)
        self.log_text.see(tk.END)
        self.log_text.config(state=tk.DISABLED)
        
        # Update the last output time
        self.last_output_time = time.time()
        
        # Filter out unnecessary or repetitive information from logs for cleaner display
        if "DEBUG:" in text or "INPUT:" in text:
            return
    
    def update_help(self, text):
        """Update help panel with guidance text."""
        self.help_text.config(state=tk.NORMAL)
        self.help_text.delete(1.0, tk.END)
        self.help_text.insert(tk.END, text)
        self.help_text.config(state=tk.DISABLED)
    
    def show_initial_help(self):
        """Show initial help message."""
        help_text = (
            "GETTING STARTED\n\n"
            "1. Enter your email address\n\n"
            "2. Enter your NCBI API key\n"
            "   (Get one at: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)\n\n"
            "3. Enter TaxID number(s) separated by commas\n"
            "   (Find TaxIDs at: https://www.ncbi.nlm.nih.gov/taxonomy)\n\n"
            "4. Click 'Save Configuration' when done.\n"
        )
        self.update_help(help_text)
    
    def show_module_1_help(self):
        """Show help for Module 1."""
        help_text = (
            "RUNNING MODULE 1\n\n"
            "Module 1 (Genetic Sequence Retrieval) downloads genetic sequences from NCBI for your species.\n\n"
            "This module performs a survey of NCBI GenBank resources using gene names and groups sequences accordingly. " 
            "First, it identifies all gene features available for your species based on the TaxID. "
            "Then it estimates how many sequences are available for each gene and presents a ranked list.\n\n" 
            "You will then be asked to select which genes you want to download by entering their numbers.\n"
            "1. Click the '1. Genetic Sequence Retrieval' button to start.\n"
            "2. The program will search for genes in your species.\n"
            "3. When prompted, select genes by entering their numbers separated by commas.\n"
            "4. For each gene, you'll be shown a list of reference sequences to choose from.\n"
            "5. After selecting a reference sequence, the program will use BLAST to find similar sequences.\n"
            "6. Only homologous regions from the BLAST search will be saved for efficient analysis."
        )
        self.update_help(help_text)
    
    def show_input_selection_help(self):
        """Show help for input selection."""
        help_text = (
            "SELECT GENES TO DOWNLOAD\n\n"
            "Option 1: Select specific genes by number\n"
            "- Enter gene numbers separated by commas\n"
            "- Use ranges like '5-7' for multiple genes\n"
            "- Example: '1,3,5-7'\n\n"
            "Option 2: Select genes by minimum sequence count\n"
            "- Enter a minimum number of sequences\n"
            "- Example: '8'\n\n"
            "Option 3: Skip downloading sequences for this taxon\n"
            "- Just enter '3' to skip"
        )
        self.update_help(help_text)
        
    def show_gene_selection_help(self):
        """Show detailed help for gene number selection."""
        help_text = (
            "SELECTING GENES TO DOWNLOAD\n\n"
            "Enter gene numbers separated by commas. You can specify:\n\n"
            "1. Individual numbers, e.g.: '1,3,8'\n"
            "   This will download genes numbered 1, 3, and 8 from the list\n\n"
            "2. Number ranges, e.g.: '5-7'\n"
            "   This will download genes numbered 5, 6, and 7\n\n"
            "3. Combinations of the above, e.g.: '1,3,5-7,10'\n"
            "   This will download genes numbered 1, 3, 5, 6, 7, and 10\n\n"
            "Example: by entering '1,3,5-7':\n"
            "- You'll select the first gene from the list\n"
            "- You'll select the third gene from the list\n"
            "- You'll select genes from fifth to seventh from the list\n\n"
            "Tip: Choose genes with higher sequence counts\n"
            "(numbers displayed next to gene names)\n"
            "for better results in subsequent modules."
        )
        self.update_help(help_text)
    
    def show_msa_groups_help(self):
        """Show help for MSA group selection."""
        help_text = (
            "SELECT SEQUENCE FILES FOR ALIGNMENT\n\n"
            "Format: numbers separated by commas and semicolons\n\n"
            "- Numbers within the same group (separated by commas) will be aligned together\n"
            "- Semicolons separate independent alignment groups\n\n"
            "Examples:\n"
            "• '0,1' - Files #0 and #1 will be aligned together\n"
            "• '0;1' - Files #0 and #1 will be aligned separately\n"
            "• '0,1,2;3,4' - Files #0, #1, #2 will be aligned together, and files #3 and #4 will form a second, independent alignment\n\n"
            "Each group will create a separate MAFFT and consensus output."
        )
        self.update_help(help_text)

    def show_reference_sequence_help(self):
        """Show help for reference sequence selection."""
        help_text = (
            "SELECT REFERENCE SEQUENCE\n\n"
            "The reference sequence serves as the query for the BLAST search to find similar sequences.\n\n"
            "Recommendations:\n"
            "• Choose a sequence of average length compared to others\n"
            "• Select a sequence from a well-documented strain/isolate when possible\n"
            "• Avoid selecting incomplete or partial sequences\n"
            "• Look for sequences with minimal ambiguous bases (N's)\n\n"
            "The selected reference sequence will be used as the BLAST query to find homologous sequences, and only matching regions will be saved for further analysis."
        )
        self.update_help(help_text)

    def show_consensus_threshold_help(self):
        """Show help for consensus threshold selection."""
        help_text = (
            "SET CONSENSUS THRESHOLD\n\n"
            "This threshold determines when a nucleotide is included in the consensus sequence:\n\n"
            "• If a nucleotide appears at or above the threshold percentage, it's included\n"
            "• If below threshold, 'N' (ambiguous) is used instead\n\n"
            "Example with 80% threshold:\n"
            "At a specific position: A A A A T (80% A, 20% T)\n"
            "Result: A (since A appears in 80% of sequences)\n\n"
            "Example with 90% threshold:\n"
            "Same position: A A A A T (80% A, 20% T)\n"
            "Result: N (since A appears in only 80% of sequences)"
        )
        self.update_help(help_text)
        
    def show_primer_amplicon_help(self):
        """Show help for primer amplicon length selection."""
        help_text = (
            "AMPLICON LENGTH RANGE\n\n"
            "Specify the desired size range for your PCR amplicons.\n\n"
            "Format: min-max (e.g., 100-250)\n\n"
            "• Shorter amplicons (70-200 bp): Better for qPCR, degraded samples\n"
            "• Medium amplicons (200-500 bp): Standard for most applications\n"
            "• Longer amplicons (500-1000 bp): For specific applications\n\n"
            "The chosen range affects amplification efficiency, specificity,\n"
            "and compatibility with your sample types."
        )
        self.update_help(help_text)
        
    def show_primer_settings_help(self):
        """Show help for primer design settings modification."""
        help_text = (
            "PRIMER DESIGN SETTINGS\n\n"
            "You are being asked whether to modify the PCR primer design parameters stored in PCR_primer_settings.txt.\n\n"
            "• If you select 'yes': The program will display current settings and pause while you edit the file manually in a text editor. You'll need to confirm when you've finished editing.\n\n"
            "• If you select 'no': The program will continue using the current settings without interruption.\n\n"
            "Modifying these settings allows you to fine-tune primer characteristics such as melting temperature, GC content, and structural constraints for your specific application."
        )
        self.update_help(help_text)

    def show_primer_sets_help(self):
        """Show help for number of primer sets to generate."""
        help_text = (
            "NUMBER OF PRIMER SETS\n\n"
            "Specify how many primer sets you want to generate.\n\n"
            "Each set includes a forward primer, reverse primer, and probe.\n\n"
            "• Recommended: 20-100 sets for most applications\n"
            "• For difficult targets: 100-500 sets may be needed\n"
            "• Maximum supported: 1000 sets\n\n"
            "Note: Generating more sets provides more options but will\n"
            "increase processing time in subsequent specificity test modules."
        )
        self.update_help(help_text)
    
    def show_specificity_settings_help(self):
        """Show help for specificity settings."""
        help_text = (
            "PRIMER SPECIFICITY SETTINGS\n\n"
            "These settings control how primers are tested against GenBank sequences:\n\n"
            "• MAX_TOTAL_MISMATCHES: Maximum mismatches allowed for a primer to match\n"
            "• MAX_3PRIME_MISMATCHES: Maximum mismatches allowed at 3' end\n"
            "• BASES_FROM_3PRIME: Size of the critical 3' region\n"
            "• MAX_MISMATCHES_PROBE: Maximum mismatches allowed for probes\n\n"
            "Lower values result in more stringent specificity tests.\n"
            "You can modify these in primer_specificity_settings.txt"
        )
        self.update_help(help_text)
    
    def show_target_taxonomy_help(self):
        """Show help for target taxonomy IDs."""
        help_text = (
            "TARGET TAXONOMY IDs\n\n"
            "NCBI Taxonomy IDs identify the species or taxonomic groups\n"
            "that you want to evaluate your primers against.\n\n"
            "• Enter one or more TaxIDs separated by commas\n"
            "• Example: '4565,4577' for specific wheat species\n"
            "• Leave empty to evaluate without taxonomic filtering\n\n"
            "TaxIDs can be found on the NCBI Taxonomy website:\n"
            "https://www.ncbi.nlm.nih.gov/taxonomy"
        )
        self.update_help(help_text)
        
    def show_module_help(self, module_num):
        """Show help for specific module."""
        module_descriptions = {
            1: "Module 1: Genetic Sequence Retrieval\n\nRetrieves genetic sequences from NCBI for your specified taxonomic IDs (TaxIDs).",
            2: "Module 2: Multiple Sequence Alignment\n\nAligns the genetic sequences using MAFFT to identify similarities and differences between sequences and generates consensus sequences.\n\nMAFFT (Multiple Alignment using Fast Fourier Transform) is a high-performance, accurate alignment tool developed by Kazutaka Katoh. MAFFT was originally published by Katoh et al. in 2002, with major updates enhancing scalability and accuracy in subsequent versions. Reference: Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution. 2013;30(4):772–780. It uses progressive and iterative refinement approaches with FFT-based scoring for rapid identification of homologous regions. MAFFT was chosen for its excellent balance of speed and accuracy compared to other alignment algorithms. The generated alignments are then used to build consensus sequences that represent the conserved regions across your input sequences.",
            3: "Module 3: PCR Primer Design\n\nDesigns PCR primers based on the consensus sequences.\n\nThis module utilizes the Primer3 library to design complete qPCR oligonucleotide sets, which include forward primer, reverse primer, and a hybridization probe for each amplicon. The design process evaluates thermodynamic parameters (Tm, GC content), secondary structure formation potential, and sequence specificity, ensuring each oligonucleotide adheres to optimal qPCR requirements. The software applies sophisticated algorithms to identify conserved regions in the consensus sequence that are suitable for species-specific amplification, then ranks primer sets based on their quality scores. Each primer set is designed to function optimally under standard qPCR conditions and is saved in a format ready for experimental validation.",
            4: "Module 4: Primer Specificity Testing\n\nTests the specificity of the designed primers against GenBank sequences to ensure they only amplify the intended targets."
        }
        
        help_text = module_descriptions.get(module_num, "No help available for this module.")
        help_text += "\n\nWhen prompted for input, please read the options carefully and enter your choice in the console input field at the bottom."
        
        self.update_help(help_text)
    
    def enable_input_controls(self):
        """Enable the input controls."""
        self.input_entry.config(state=tk.NORMAL)
        self.submit_button.config(state=tk.NORMAL)
        self.input_entry.focus_set()
        
    def disable_input_controls(self):
        """Disable the input controls."""
        self.input_entry.config(state=tk.DISABLED)
        self.submit_button.config(state=tk.DISABLED)
    
    def submit_input(self):
        """Submit user input using file-based approach."""
        if not self.running or not self.waiting_for_input:
            return
            
        user_input = self.input_var.get()
        if not user_input:
            return
            
        # Write the response to the response file
        input_response_file = os.path.join(os.getcwd(), 'input_response.txt')
        
        try:
            with open(input_response_file, 'w', encoding='utf-8') as f:
                f.write(user_input)
                
            # Log input in the output window
            self.update_log(f"> {user_input}\n")
            
            # Clear input field and hide the message
            self.input_var.set("")
            self.input_message.config(text="")
            self.input_message.pack_forget()
            
            # Disable input until next request
            self.waiting_for_input = False
            self.disable_input_controls()
            
            self.status_var.set("Running...")
        except Exception as e:
            self.update_log(f"Error writing input response: {e}\n")
    
    # This function is no longer needed with file-based input
    # but keeping it as a stub for compatibility
    def input_writer(self, process):
        """Thread function for backwards compatibility."""
        try:
            while process.poll() is None:
                time.sleep(0.5)
        except Exception as e:
            self.output_queue.put(f"Input thread error: {e}\n")
        finally:
            if process.stdin:
                process.stdin.close()
    
    def clear_log(self):
        """Clear the log field."""
        self.log_text.config(state=tk.NORMAL)
        self.log_text.delete(1.0, tk.END)
        self.log_text.config(state=tk.DISABLED)
    
    def read_output(self):
        """Read output from queue and update logs."""
        if self.running:
            try:
                while True:
                    line = self.output_queue.get_nowait()
                    
                    # Filter out redundant or noisy log entries
                    if "records processed" in line and "in this batch" in line:
                        # Skip redundant processing messages
                        pass
                    elif "Processing batch" in line and "overall progress" in line:
                        # Keep progress updates
                        self.update_log(line)
                        # Update status bar with progress information
                        progress_match = re.search(r'progress: (\d+\.\d+)%', line)
                        if progress_match:
                            progress = float(progress_match.group(1))
                            self.status_var.set(f"Running - {progress:.1f}% complete")
                    elif any(x in line for x in [
                        "Fetching and analyzing records", 
                        "DEBUG:", 
                        "successfully saved", 
                        "retrying"]):
                        # Skip verbose messages
                        pass
                    else:
                        # Keep all other messages
                        self.update_log(line)
                    
                    self.output_queue.task_done()
            except queue.Empty:
                pass
            
            # Check for input request file
            self.check_for_input_requests()
            
            # Check if process is still running
            if self.process and self.process.poll() is not None:
                returncode = self.process.poll()
                if returncode == 0:
                    self.update_log("\n--- Process completed successfully ---\n")
                else:
                    self.update_log(f"\n--- Process ended with error (code: {returncode}) ---\n")
                
                # Cancel timeout timer
                if self.timeout_id:
                    self.root.after_cancel(self.timeout_id)
                    self.timeout_id = None
                
                # Cancel input checker
                if self.input_checker_id:
                    self.root.after_cancel(self.input_checker_id)
                    self.input_checker_id = None
                
                # Remove any input request/response files
                self.cleanup_input_files()
                
                # Hide the input message if it's visible
                self.input_message.pack_forget()
                
                self.running = False
                self.waiting_for_input = False
                self.disable_input_controls()
                self.status_var.set("Ready")
                self.progress_bar.stop()
                self.enable_buttons()
                
                # Show help for next module
                next_module = self.current_module + 1
                if 1 <= next_module <= 8:
                    self.show_module_help(next_module)
            else:
                # Program is still running
                self.root.after(100, self.read_output)
    
    def check_for_input_requests(self):
        """Check if the script is requesting input via file"""
        input_request_file = os.path.join(os.getcwd(), 'input_request.txt')
        
        if os.path.exists(input_request_file):
            try:
                with open(input_request_file, 'r', encoding='utf-8') as f:
                    prompt = f.read().strip()
                
                # Enable input and set message with yellow highlight
                self.waiting_for_input = True
                self.enable_input_controls()
                
                # Show the input message with yellow background
                self.input_message.config(text=f"Input requested: {prompt}")
                self.input_message.pack(fill=tk.X, pady=2)
                
                self.status_var.set("Waiting for input. Please provide a response.")
                
                # Show appropriate help based on the prompt
                if self.current_module == 1 and "Enter gene numbers separated by commas" in prompt:
                    self.show_gene_selection_help()
                elif self.current_module == 1 and "Please select one of the following options" in prompt:
                    self.show_input_selection_help()
                elif self.current_module == 2 and "Enter groups of file indices" in prompt:
                    self.show_msa_groups_help()
                elif self.current_module == 2 and "Choose reference sequence" in prompt:
                    self.show_reference_sequence_help()
                elif self.current_module == 2 and "Enter consensus threshold percentage" in prompt:
                    self.show_consensus_threshold_help()
                elif self.current_module == 3 and "Do you want to modify primer design settings" in prompt:
                    self.show_primer_settings_help()
                elif self.current_module == 3 and "Enter amplicon length range" in prompt:
                    self.show_primer_amplicon_help()
                elif self.current_module == 3 and "Enter maximum number of primer sets" in prompt:
                    self.show_primer_sets_help()
                elif self.current_module == 4 and "modify specificity test settings" in prompt.lower():
                    self.show_specificity_settings_help()
                elif self.current_module == 4 and "Enter NCBI Taxonomy ID" in prompt:
                    self.show_target_taxonomy_help()
            except Exception as e:
                self.update_log(f"Error reading input request: {e}\n")
    
    def cleanup_input_files(self):
        """Remove any input request/response files"""
        input_request_file = os.path.join(os.getcwd(), 'input_request.txt')
        input_response_file = os.path.join(os.getcwd(), 'input_response.txt')
        
        for file in [input_request_file, input_response_file]:
            if os.path.exists(file):
                try:
                    os.remove(file)
                except Exception as e:
                    self.update_log(f"Error removing file {file}: {e}\n")
    
    def output_reader(self, pipe, queue):
        """Function to read process output and place it in queue."""
        try:
            for line in iter(pipe.readline, b''):
                text = line.decode('utf-8', errors='replace')
                queue.put(text)
                time.sleep(0.01)  # Small delay to avoid overloading the queue
        except Exception as e:
            queue.put(f"Output reading error: {e}\n")
        finally:
            pipe.close()

    def disable_buttons(self):
        """Disable all module buttons during processing."""
        for button in self.module_buttons:
            button.config(state=tk.DISABLED)
    
    def enable_buttons(self):
        """Enable all module buttons after processing completes."""
        for button in self.module_buttons:
            button.config(state=tk.NORMAL)

    def check_requirements(self, module_num):
        """Check if requirements for the module are met."""
        # Check NCBI data
        if module_num >= 1 and (not self.email.get() or not self.api_key.get()):
            messagebox.showerror("Configuration Error", "You must provide email and API key.")
            return False
            
        # Check directories for previous modules
        for i in range(1, module_num):
            directory = Path(f"{i}_")
            if not directory.exists() or not any(directory.iterdir()):
                if not messagebox.askyesno("Missing Data", 
                                          f"Directory '{directory}' does not exist or is empty. "
                                          f"Module {i} was not completed. Do you want to continue anyway?"):
                    return False
        
        # Check TaxID for module 8
        if module_num == 8 and not self.taxids.get():
            messagebox.showerror("Missing TaxID", "You must provide TaxID numbers for Module 8.")
            return False
            
        return True
    
    def save_user_data(self):
        """Save user data to files."""
        os.makedirs('user_data', exist_ok=True)
        
        try:
            with open(os.path.join('user_data', 'email.txt'), 'w') as f:
                f.write(self.email.get())
                
            with open(os.path.join('user_data', 'api_key.txt'), 'w') as f:
                f.write(self.api_key.get())
                
            return True
        except Exception as e:
            self.update_log(f"Error saving user data: {e}\n")
            return False
    
    def get_module_script(self, module_num):
        """Return script name for the given module."""
        module_scripts = {
            1: "1_Genetic_Background.py",
            2: "2_MSA_Alignment.py",
            3: "3_PCR_Primers_Design.py",
            4: "4_Primers_Specificity.py"
        }
        return module_scripts.get(module_num)

    def run_module(self, module_num):
        """Run selected program module."""
        if self.running:
            messagebox.showwarning("Process Running", "Another module is already running. Stop it before starting a new one.")
            return
            
        if not self.check_requirements(module_num):
            return
            
        script_name = self.get_module_script(module_num)
        if not script_name or not os.path.exists(script_name):
            messagebox.showerror("Error", f"Script for module {module_num} not found")
            return
            
        # Save user data to files
        if not self.save_user_data():
            if not messagebox.askyesno("Error", "Failed to save user data. Continue anyway?"):
                return
                
        self.save_config()
        
        # Create required directories
        os.makedirs(f"{module_num}_", exist_ok=True)
        
        # Clean up any stray input files
        self.cleanup_input_files()
        
        # Run the module
        self.current_module = module_num
        self.running = False  # Set to True after process starts
        self.waiting_for_input = False
        
        # Show help for current module
        self.show_module_help(module_num)
        
        # Clear the log
        self.clear_log()
        
        # Disable buttons while running
        self.disable_buttons()
        
        # Start progress bar
        self.progress_bar.start(10)
        
        self.update_log(f"Starting module {module_num}: {script_name}\n")
        self.update_log("=" * 50 + "\n")
        
        self.status_var.set(f"Running module {module_num}")
        
        # Record start time
        self.last_output_time = time.time()
        
        # Create a timestamp for log files
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        log_dir = "logs"
        os.makedirs(log_dir, exist_ok=True)
        log_file = os.path.join(log_dir, f"module{module_num}_{timestamp}.log")
        
        try:
            # Modified Python environment - ensure stdin isn't buffered 
            env = os.environ.copy()
            env['PYTHONUNBUFFERED'] = '1'
            
            # Run in interactive mode to enable communication
            cmd = [sys.executable, script_name, "--interactive"]
            
            # Set up pipes for process communication
            self.process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                stdin=subprocess.PIPE,
                bufsize=1,  # Line buffered
                universal_newlines=False,
                env=env
            )
            
            # Start threads for reading output
            t_output = threading.Thread(target=self.output_reader, args=(self.process.stdout, self.output_queue))
            t_output.daemon = True
            t_output.start()
            
            # Start thread for writing input (for backward compatibility)
            t_input = threading.Thread(target=self.input_writer, args=(self.process,))
            t_input.daemon = True
            t_input.start()
            
            # Start checking for input requests
            self.input_checker_id = self.root.after(500, self.check_for_input_requests)
            
            # Disable input controls initially
            self.disable_input_controls()
            
            # Start reading output
            self.running = True
            self.read_output()
            
        except Exception as e:
            self.update_log(f"Error starting process: {e}\n")
            self.running = False
            self.waiting_for_input = False
            self.status_var.set("Error")
            self.progress_bar.stop()
            self.enable_buttons()
    
    def stop_process(self):
        """Stop currently running process."""
        if not self.running or not self.process:
            messagebox.showinfo("Information", "No process is currently running.")
            return
        
        try:
            # Try to terminate process gracefully first
            if sys.platform == 'win32':
                self.process.terminate()
            else:
                # On Unix systems, try SIGTERM first
                self.process.send_signal(signal.SIGTERM)
            
            # Give process time to close
            for i in range(5):  # Check for 5 seconds
                if self.process.poll() is not None:
                    break
                self.update_log(".")
                time.sleep(1)
                
            # If process is still running, force kill it
            if self.process.poll() is None:
                if sys.platform == 'win32':
                    self.process.kill()
                else:
                    self.process.send_signal(signal.SIGKILL)
                self.update_log("\nProcess had to be forcefully terminated.\n")
            else:
                self.update_log("\nProcess terminated successfully.\n")
            
            # Cancel timeout timer
            if self.timeout_id:
                self.root.after_cancel(self.timeout_id)
                self.timeout_id = None
                
            # Cancel input checker
            if self.input_checker_id:
                self.root.after_cancel(self.input_checker_id)
                self.input_checker_id = None
                
            # Clean up any input files
            self.cleanup_input_files()
            
            self.running = False
            self.waiting_for_input = False
            self.disable_input_controls()
            self.input_message.config(text="")
            self.input_message.pack_forget()
            self.status_var.set("Process stopped")
            self.progress_bar.stop()
            self.enable_buttons()
            
        except Exception as e:
            self.update_log(f"\nError stopping process: {e}\n")
    
    def restart_application(self):
        """Restart the entire application."""
        if self.running:
            if not messagebox.askyesno("Warning", "A process is currently running. Stop it and restart the application?"):
                return
            self.stop_process()
        
        self.update_log("\nRestarting application...\n")
        
        # Clean up any input files before restart
        self.cleanup_input_files()
        
        # Python executable and this script path
        python = sys.executable
        script = os.path.abspath(sys.argv[0])
        
        # Close the current application
        self.root.destroy()
        
        # Start a new instance
        os.execl(python, python, script)


def main():
    root = tk.Tk()
    
    # Set application icon if available
    try:
        icon_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "icon.ico")
        if os.path.exists(icon_path):
            root.iconbitmap(icon_path)
    except:
        pass
        
    # Configure styles for input message with yellow background
    style = ttk.Style()
    style.configure("YellowBg.TLabel", background="#FFFF99", foreground="black")
    
    app = PrimeSpecPCRApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
