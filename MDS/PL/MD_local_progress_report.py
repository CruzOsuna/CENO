# MDS Progress Report
import os
import logging
import matplotlib.pyplot as plt
from datetime import datetime

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Function to detect all folders ending with "_MDS"
def get_protein_folders():
    return [folder for folder in os.listdir() if os.path.isdir(folder) and folder.endswith("_MDS")]

# Function to check if the Energy Minimization (EM) has been successfully completed
def is_em_completed(folder):
    em_log_file = os.path.join(folder, "EM.log")
    if os.path.exists(em_log_file):
        with open(em_log_file, 'r') as log_file:
            log_content = log_file.read()
            return "Finished mdrun on rank 0" in log_content
    return False

# Function to check if the NVT equilibration has been successfully completed
def is_nvt_equilibration_completed(folder):
    nvt_log_file = os.path.join(folder, "NVT.log")
    if os.path.exists(nvt_log_file):
        with open(nvt_log_file, 'r') as log_file:
            log_content = log_file.read()
            return "Finished mdrun on rank 0" in log_content
    return False

# Function to check if the NPT equilibration has been successfully completed
def is_npt_equilibration_completed(folder):
    npt_log_file = os.path.join(folder, "NPT.log")
    if os.path.exists(npt_log_file):
        with open(npt_log_file, 'r') as log_file:
            log_content = log_file.read()
            return "Finished mdrun on rank 0" in log_content
    return False

# Function to check if the MD simulation has been successfully completed
def is_md_simulation_completed(folder):
    md_log_file = os.path.join(folder, "analisis", "MD.log")
    if os.path.exists(md_log_file):
        with open(md_log_file, 'r') as log_file:
            log_content = log_file.read()
            return "Finished mdrun on rank 0" in log_content
    return False

# Function to calculate the size of each protein folder
def get_folder_size(folder):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(folder):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size / (1024 ** 2)  # Convert to MB

# Function to generate the monitoring report and plots
def generate_monitoring_report_and_plots():
    protein_folders = get_protein_folders()
    total_proteins = len(protein_folders)
    
    em_completed = []
    nvt_completed = []
    npt_completed = []
    md_completed = []
    protein_progress = {}
    folder_sizes = {}

    for folder in protein_folders:
        em_done = is_em_completed(folder)
        nvt_done = is_nvt_equilibration_completed(folder)
        npt_done = is_npt_equilibration_completed(folder)
        md_done = is_md_simulation_completed(folder)

        if em_done:
            em_completed.append(folder)
        if nvt_done:
            nvt_completed.append(folder)
        if npt_done:
            npt_completed.append(folder)
        if md_done:
            md_completed.append(folder)

        # Calculate the total progress for each protein
        progress = (1/12 if em_done else 0) + (1/12 if nvt_done else 0) + (1/12 if npt_done else 0) + (9/12 if md_done else 0)
        protein_progress[folder] = progress * 100

        # Get the size of each protein folder
        folder_sizes[folder] = get_folder_size(folder)

    # Calculate global percentages
    em_percentage = (len(em_completed) / total_proteins) * 100 if total_proteins > 0 else 0
    nvt_percentage = (len(nvt_completed) / total_proteins) * 100 if total_proteins > 0 else 0
    npt_percentage = (len(npt_completed) / total_proteins) * 100 if total_proteins > 0 else 0
    md_percentage = (len(md_completed) / total_proteins) * 100 if total_proteins > 0 else 0

    # Plot 1: Global progress percentage for each step
    plt.figure(figsize=(12, 8))
    steps = ['EM', 'NVT', 'NPT', 'MD']
    percentages = [em_percentage, nvt_percentage, npt_percentage, md_percentage]
    plt.bar(steps, percentages)
    plt.title('Global Progress Percentage per Step')
    plt.ylabel('Percentage (%)')
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig('global_progress.png')
    plt.close()

    # Plot 2: Progress percentage for each individual protein
    plt.figure(figsize=(12, 8))  # Increase the figure size for better readability
    plt.barh(list(protein_progress.keys()), list(protein_progress.values()))
    plt.title('Progress Percentage per Protein')
    plt.xlabel('Percentage (%)')
    plt.xticks(fontsize=10)  # Adjust fontsize for x-axis
    plt.yticks(fontsize=8)  # Adjust fontsize for y-axis
    plt.xlim(0, 100)
    plt.tight_layout()
    plt.savefig('progress_per_protein.png')
    plt.close()

    # Plot 3: Storage size of each protein folder
    plt.figure(figsize=(12, 8))  # Increase the figure size for better readability
    plt.barh(list(folder_sizes.keys()), list(folder_sizes.values()))
    plt.title('Storage Size per Protein')
    plt.xlabel('Size (MB)')
    plt.xticks(fontsize=10)  # Adjust fontsize for x-axis
    plt.yticks(fontsize=8)  # Adjust fontsize for y-axis
    plt.tight_layout()
    plt.savefig('storage_per_protein.png')
    plt.close()

    # Generate the monitoring report
    report_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    report_content = (
        f"Monitoring Report - {report_time}\n"
        f"Total proteins: {total_proteins}\n\n"
        f"EM Progress: {len(em_completed)} completed ({em_percentage:.2f}%)\n"
        f"NVT Equilibration Progress: {len(nvt_completed)} completed ({nvt_percentage:.2f}%)\n"
        f"NPT Equilibration Progress: {len(npt_completed)} completed ({npt_percentage:.2f}%)\n"
        f"MD Simulation Progress: {len(md_completed)} completed ({md_percentage:.2f}%)\n\n"
        f"Proteins that completed EM:\n" + "\n".join(em_completed) + "\n\n"
        f"Proteins that completed NVT:\n" + "\n".join(nvt_completed) + "\n\n"
        f"Proteins that completed NPT:\n" + "\n".join(npt_completed) + "\n\n"
        f"Proteins that completed MD:\n" + "\n".join(md_completed) + "\n\n"
    )

    # Save the report to a file
    report_filename = "complete_monitoring_report.txt"
    with open(report_filename, "w") as report_file:
        report_file.write(report_content)
    
    logging.info(f"Monitoring report and plots saved. Report: {report_filename}")

# Run the monitoring and plotting function
if __name__ == "__main__":
    generate_monitoring_report_and_plots()
