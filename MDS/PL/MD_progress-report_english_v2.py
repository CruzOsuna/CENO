#MD_progress-report_v2

#English version

import logging
from pathlib import Path
import matplotlib.pyplot as plt
from datetime import datetime
import concurrent.futures
import smtplib
import ssl
from email.message import EmailMessage
import mimetypes
import schedule
import time
import imaplib
import email
import getpass
import argparse

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("monitoring.log"),
        logging.StreamHandler()
    ]
)

def get_protein_folders(input_dir):
    """
    Detects all folders in the specified directory that end with '_MDS'.

    Args:
        input_dir (Path): The directory where protein folders are located.

    Returns:
        list of Path: List of Path objects corresponding to protein folders.
    """
    return [folder for folder in input_dir.iterdir() if folder.is_dir() and folder.name.endswith("_MDS")]

def is_step_completed(folder, log_filename):
    """
    Checks if a simulation step has been successfully completed by looking for a specific message in the log file.

    Args:
        folder (Path): The protein folder.
        log_filename (str): The name of the log file for the simulation step.

    Returns:
        bool: True if the step is completed, False otherwise.
    """
    log_file = folder / log_filename
    if log_file.exists():
        try:
            with log_file.open('r', encoding='utf-8') as f:
                log_content = f.read()
                return "Finished mdrun on rank 0" in log_content
        except Exception as e:
            logging.error(f"Error reading {log_file}: {e}")
            return False
    return False

def get_folder_size(folder):
    """
    Calculates the size of a folder in MB.

    Args:
        folder (Path): The folder whose size is to be calculated.

    Returns:
        float: The size of the folder in MB.
    """
    total_size = 0
    try:
        for f in folder.rglob('*'):
            if f.is_file():
                total_size += f.stat().st_size
    except Exception as e:
        logging.error(f"Error calculating the size for {folder}: {e}")
    return total_size / (1024 ** 2)  # Convert to MB

def process_folder(folder):
    """
    Processes a protein folder to determine the completion status of each simulation step and calculate the folder size.

    Args:
        folder (Path): The protein folder.

    Returns:
        tuple: Contains the folder name, completion status of steps, progress, and folder size.
    """
    em_done = is_step_completed(folder, "EM.log")
    nvt_done = is_step_completed(folder, "NVT.log")
    npt_done = is_step_completed(folder, "NPT.log")
    md_done = is_step_completed(folder / "analisis", "MD.log")

    # Calculate progress
    progress = (1/12 if em_done else 0) + \
               (1/12 if nvt_done else 0) + \
               (1/12 if npt_done else 0) + \
               (9/12 if md_done else 0)
    progress *= 100  # Convert to percentage

    # Get folder size
    folder_size = get_folder_size(folder)

    return (folder.name, em_done, nvt_done, npt_done, md_done, progress, folder_size)

def generate_monitoring_report_and_plots(input_dir, output_dir):
    """
    Generates the monitoring report and plots for molecular dynamics simulations.

    Args:
        input_dir (Path): The input directory containing protein folders.
        output_dir (Path): The output directory for saving reports and plots.
    """
    protein_folders = get_protein_folders(input_dir)
    total_proteins = len(protein_folders)

    if total_proteins == 0:
        logging.warning("No protein folders found.")
        return

    # Parallel processing
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(executor.map(process_folder, protein_folders))

    # Initialize lists and dictionaries to store results
    em_completed = []
    nvt_completed = []
    npt_completed = []
    md_completed = []
    protein_progress = {}
    folder_sizes = {}

    # Process results
    for result in results:
        folder_name, em_done, nvt_done, npt_done, md_done, progress, folder_size = result

        if em_done:
            em_completed.append(folder_name)
        if nvt_done:
            nvt_completed.append(folder_name)
        if npt_done:
            npt_completed.append(folder_name)
        if md_done:
            md_completed.append(folder_name)

        protein_progress[folder_name] = progress
        folder_sizes[folder_name] = folder_size

    # Calculate global percentages
    em_percentage = (len(em_completed) / total_proteins) * 100
    nvt_percentage = (len(nvt_completed) / total_proteins) * 100
    npt_percentage = (len(npt_completed) / total_proteins) * 100
    md_percentage = (len(md_completed) / total_proteins) * 100

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Plot 1: Global progress percentage for each step
    plt.figure(figsize=(12, 8))
    steps = ['EM', 'NVT', 'NPT', 'MD']
    percentages = [em_percentage, nvt_percentage, npt_percentage, md_percentage]
    plt.bar(steps, percentages)
    plt.title('Global Progress Percentage per Step')
    plt.ylabel('Percentage (%)')
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(output_dir / 'global_progress.png')
    plt.close()

    # Plot 2: Progress percentage for each individual protein
    plt.figure(figsize=(12, 8))
    # Sort proteins by progress
    sorted_proteins = sorted(protein_progress.items(), key=lambda x: x[1])
    protein_names = [x[0] for x in sorted_proteins]
    progress_values = [x[1] for x in sorted_proteins]
    plt.barh(protein_names, progress_values)
    plt.title('Progress Percentage per Protein')
    plt.xlabel('Percentage (%)')
    plt.xlim(0, 100)
    plt.tight_layout()
    plt.savefig(output_dir / 'progress_per_protein.png')
    plt.close()

    # Plot 3: Storage size of each protein folder
    plt.figure(figsize=(12, 8))
    # Sort proteins by folder size
    sorted_sizes = sorted(folder_sizes.items(), key=lambda x: x[1])
    protein_names_size = [x[0] for x in sorted_sizes]
    size_values = [x[1] for x in sorted_sizes]
    plt.barh(protein_names_size, size_values)
    plt.title('Storage Size per Protein')
    plt.xlabel('Size (MB)')
    plt.tight_layout()
    plt.savefig(output_dir / 'storage_per_protein.png')
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
    report_filename = output_dir / "complete_monitoring_report.txt"
    with report_filename.open("w", encoding='utf-8') as report_file:
        report_file.write(report_content)

    logging.info(f"Monitoring report and plots saved in {output_dir}")

def send_email_report(sender_email, sender_password, recipient_email, subject, body, attachments):
    """
    Sends an email with the specified attachments.

    Args:
        sender_email (str): Sender's email address.
        sender_password (str): Sender's email password.
        recipient_email (str): Recipient's email address.
        subject (str): Email subject.
        body (str): Email body.
        attachments (list): List of file paths to attach.
    """
    # Create the email message
    msg = EmailMessage()
    msg['From'] = sender_email
    msg['To'] = recipient_email
    msg['Subject'] = subject
    msg.set_content(body)

    # Attach files
    for file in attachments:
        try:
            mime_type, _ = mimetypes.guess_type(file)
            if mime_type is None:
                mime_type = 'application/octet-stream'
            mime_type, mime_subtype = mime_type.split('/')

            with open(file, 'rb') as f:
                msg.add_attachment(f.read(),
                                   maintype=mime_type,
                                   subtype=mime_subtype,
                                   filename=Path(file).name)
        except Exception as e:
            logging.error(f"Error attaching file {file}: {e}")

    # Send the email
    try:
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL('smtp.gmail.com', 465, context=context) as server:
            server.login(sender_email, sender_password)
            server.send_message(msg)
        logging.info("Email sent successfully.")
    except Exception as e:
        logging.error(f"Error sending email: {e}")

def check_email_for_request(email_account, email_password, search_subject):
    """
    Checks for unread emails with a specific subject to trigger report generation.

    Args:
        email_account (str): The email account to check.
        email_password (str): The password for the email account.
        search_subject (str): The subject line to look for.

    Returns:
        bool: True if the request email is found, False otherwise.
    """
    try:
        mail = imaplib.IMAP4_SSL('imap.gmail.com')
        mail.login(email_account, email_password)
        mail.select('inbox')

        # Search for unread emails with the specified subject
        result, data = mail.search(None, f'(UNSEEN SUBJECT "{search_subject}")')

        mail_ids = data[0].split()
        if mail_ids:
            # Optionally mark the emails as read
            for mail_id in mail_ids:
                mail.store(mail_id, '+FLAGS', '\\Seen')
            return True
        return False
    except Exception as e:
        logging.error(f"Error checking email for request: {e}")
        return False

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Molecular Dynamics Simulation Monitoring Script')
    parser.add_argument('--input_dir', type=str, default='.', help='Input directory containing protein folders.')
    parser.add_argument('--output_dir', type=str, default='.', help='Output directory for reports and plots.')
    parser.add_argument('--sender_email', type=str, required=True, help="Sender's email address.")
    parser.add_argument('--recipient_email', type=str, required=True, help="Recipient's email address.")
    parser.add_argument('--email_password', type=str, help="Sender's email password.")
    parser.add_argument('--interval_hours', type=int, default=6, help='Interval in hours for periodic reporting.')
    parser.add_argument('--email_subject', type=str, default='MD Report Request', help='Subject line to trigger report generation.')
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    sender_email = args.sender_email
    recipient_email = args.recipient_email
    email_password = args.email_password
    if email_password is None:
        email_password = getpass.getpass(prompt="Email password: ")
    interval_hours = args.interval_hours
    email_subject = args.email_subject

    # Function to generate the report and send the email
    def job():
        generate_monitoring_report_and_plots(input_dir, output_dir)

        # Prepare the email details
        subject = 'Molecular Dynamics Simulation Monitoring Report'
        body = 'Please find the attached monitoring report and generated plots.'
        attachments = [
            str(output_dir / 'complete_monitoring_report.txt'),
            str(output_dir / 'global_progress.png'),
            str(output_dir / 'progress_per_protein.png'),
            str(output_dir / 'storage_per_protein.png')
        ]

        send_email_report(sender_email, email_password, recipient_email, subject, body, attachments)
        logging.info(f"Monitoring report and plots sent to {recipient_email}.")

    # Schedule the periodic task
    schedule.every(interval_hours).hours.do(job)

    # Start the main loop
    while True:
        # Check for email request
        if check_email_for_request(sender_email, email_password, email_subject):
            logging.info("Report request received via email.")
            job()
        schedule.run_pending()
        time.sleep(60)  # Wait 60 seconds before checking again

if __name__ == "__main__":
    main()
