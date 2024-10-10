#MD_progress-report_v2

#Version en español 

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

# Configuración de logging
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
    Detecta todas las carpetas en el directorio especificado que terminan con '_MDS'.

    Args:
        input_dir (Path): El directorio donde buscar las carpetas de proteínas.

    Returns:
        list of Path: Lista de objetos Path correspondientes a las carpetas de proteínas.
    """
    return [folder for folder in input_dir.iterdir() if folder.is_dir() and folder.name.endswith("_MDS")]

def is_step_completed(folder, log_filename):
    """
    Verifica si un paso de simulación se ha completado exitosamente buscando un mensaje específico en el archivo de log.

    Args:
        folder (Path): La carpeta de la proteína.
        log_filename (str): El nombre del archivo de log para el paso de simulación.

    Returns:
        bool: True si el paso está completado, False en caso contrario.
    """
    log_file = folder / log_filename
    if log_file.exists():
        try:
            with log_file.open('r', encoding='utf-8') as f:
                log_content = f.read()
                return "Finished mdrun on rank 0" in log_content
        except Exception as e:
            logging.error(f"Error leyendo {log_file}: {e}")
            return False
    return False

def get_folder_size(folder):
    """
    Calcula el tamaño de una carpeta en MB.

    Args:
        folder (Path): La carpeta de la cual calcular el tamaño.

    Returns:
        float: El tamaño de la carpeta en MB.
    """
    total_size = 0
    try:
        for f in folder.rglob('*'):
            if f.is_file():
                total_size += f.stat().st_size
    except Exception as e:
        logging.error(f"Error calculando el tamaño para {folder}: {e}")
    return total_size / (1024 ** 2)  # Convertir a MB

def process_folder(folder):
    """
    Procesa una carpeta de proteína para determinar el estado de finalización de cada paso de simulación y calcular el tamaño de la carpeta.

    Args:
        folder (Path): La carpeta de la proteína.

    Returns:
        tuple: Contiene el nombre de la carpeta, estado de finalización de los pasos, progreso y tamaño de la carpeta.
    """
    em_done = is_step_completed(folder, "EM.log")
    nvt_done = is_step_completed(folder, "NVT.log")
    npt_done = is_step_completed(folder, "NPT.log")
    md_done = is_step_completed(folder / "analisis", "MD.log")

    # Calcular el progreso
    progress = (1/12 if em_done else 0) + \
               (1/12 if nvt_done else 0) + \
               (1/12 if npt_done else 0) + \
               (9/12 if md_done else 0)
    progress *= 100  # Convertir a porcentaje

    # Obtener el tamaño de la carpeta
    folder_size = get_folder_size(folder)

    return (folder.name, em_done, nvt_done, npt_done, md_done, progress, folder_size)

def generate_monitoring_report_and_plots(input_dir, output_dir):
    """
    Genera el informe de monitoreo y los gráficos para las simulaciones de dinámica molecular.

    Args:
        input_dir (Path): El directorio de entrada que contiene las carpetas de proteínas.
        output_dir (Path): El directorio de salida para guardar los informes y gráficos.
    """
    protein_folders = get_protein_folders(input_dir)
    total_proteins = len(protein_folders)

    if total_proteins == 0:
        logging.warning("No se encontraron carpetas de proteínas.")
        return

    # Procesamiento paralelo
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(executor.map(process_folder, protein_folders))

    # Inicializar listas y diccionarios para almacenar resultados
    em_completed = []
    nvt_completed = []
    npt_completed = []
    md_completed = []
    protein_progress = {}
    folder_sizes = {}

    # Procesar resultados
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

    # Calcular porcentajes globales
    em_percentage = (len(em_completed) / total_proteins) * 100
    nvt_percentage = (len(nvt_completed) / total_proteins) * 100
    npt_percentage = (len(npt_completed) / total_proteins) * 100
    md_percentage = (len(md_completed) / total_proteins) * 100

    # Crear el directorio de salida si no existe
    output_dir.mkdir(parents=True, exist_ok=True)

    # Gráfico 1: Porcentaje de progreso global para cada paso
    plt.figure(figsize=(12, 8))
    steps = ['EM', 'NVT', 'NPT', 'MD']
    percentages = [em_percentage, nvt_percentage, npt_percentage, md_percentage]
    plt.bar(steps, percentages)
    plt.title('Porcentaje de Progreso Global por Paso')
    plt.ylabel('Porcentaje (%)')
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(output_dir / 'global_progress.png')
    plt.close()

    # Gráfico 2: Porcentaje de progreso para cada proteína individual
    plt.figure(figsize=(12, 8))
    # Ordenar proteínas por progreso
    sorted_proteins = sorted(protein_progress.items(), key=lambda x: x[1])
    protein_names = [x[0] for x in sorted_proteins]
    progress_values = [x[1] for x in sorted_proteins]
    plt.barh(protein_names, progress_values)
    plt.title('Porcentaje de Progreso por Proteína')
    plt.xlabel('Porcentaje (%)')
    plt.xlim(0, 100)
    plt.tight_layout()
    plt.savefig(output_dir / 'progress_per_protein.png')
    plt.close()

    # Gráfico 3: Tamaño de almacenamiento de cada carpeta de proteína
    plt.figure(figsize=(12, 8))
    # Ordenar proteínas por tamaño de carpeta
    sorted_sizes = sorted(folder_sizes.items(), key=lambda x: x[1])
    protein_names_size = [x[0] for x in sorted_sizes]
    size_values = [x[1] for x in sorted_sizes]
    plt.barh(protein_names_size, size_values)
    plt.title('Tamaño de Almacenamiento por Proteína')
    plt.xlabel('Tamaño (MB)')
    plt.tight_layout()
    plt.savefig(output_dir / 'storage_per_protein.png')
    plt.close()

    # Generar el informe de monitoreo
    report_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    report_content = (
        f"Informe de Monitoreo - {report_time}\n"
        f"Total de proteínas: {total_proteins}\n\n"
        f"Progreso EM: {len(em_completed)} completado ({em_percentage:.2f}%)\n"
        f"Progreso de Equilibración NVT: {len(nvt_completed)} completado ({nvt_percentage:.2f}%)\n"
        f"Progreso de Equilibración NPT: {len(npt_completed)} completado ({npt_percentage:.2f}%)\n"
        f"Progreso de Simulación MD: {len(md_completed)} completado ({md_percentage:.2f}%)\n\n"
        f"Proteínas que completaron EM:\n" + "\n".join(em_completed) + "\n\n"
        f"Proteínas que completaron NVT:\n" + "\n".join(nvt_completed) + "\n\n"
        f"Proteínas que completaron NPT:\n" + "\n".join(npt_completed) + "\n\n"
        f"Proteínas que completaron MD:\n" + "\n".join(md_completed) + "\n\n"
    )

    # Guardar el informe en un archivo
    report_filename = output_dir / "complete_monitoring_report.txt"
    with report_filename.open("w", encoding='utf-8') as report_file:
        report_file.write(report_content)

    logging.info(f"Informe de monitoreo y gráficos guardados en {output_dir}")

def send_email_report(sender_email, sender_password, recipient_email, subject, body, attachments):
    """
    Envía un correo electrónico con los archivos adjuntos especificados.

    Args:
        sender_email (str): Correo electrónico del remitente.
        sender_password (str): Contraseña del correo electrónico del remitente.
        recipient_email (str): Correo electrónico del destinatario.
        subject (str): Asunto del correo electrónico.
        body (str): Cuerpo del correo electrónico.
        attachments (list): Lista de rutas a archivos que se adjuntarán.
    """
    # Crear el mensaje de correo electrónico
    msg = EmailMessage()
    msg['From'] = sender_email
    msg['To'] = recipient_email
    msg['Subject'] = subject
    msg.set_content(body)

    # Adjuntar archivos
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
            logging.error(f"Error adjuntando el archivo {file}: {e}")

    # Enviar el correo electrónico
    try:
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL('smtp.gmail.com', 465, context=context) as server:
            server.login(sender_email, sender_password)
            server.send_message(msg)
        logging.info("Correo electrónico enviado exitosamente.")
    except Exception as e:
        logging.error(f"Error al enviar el correo electrónico: {e}")

def check_email_for_request(email_account, email_password, search_subject):
    """
    Verifica si hay correos electrónicos no leídos con un asunto específico para desencadenar la generación del informe.

    Args:
        email_account (str): La cuenta de correo electrónico a verificar.
        email_password (str): La contraseña para la cuenta de correo electrónico.
        search_subject (str): La línea de asunto a buscar.

    Returns:
        bool: True si se encuentra el correo de solicitud, False en caso contrario.
    """
    try:
        mail = imaplib.IMAP4_SSL('imap.gmail.com')
        mail.login(email_account, email_password)
        mail.select('inbox')

        # Buscar correos no leídos con el asunto especificado
        result, data = mail.search(None, f'(UNSEEN SUBJECT "{search_subject}")')

        mail_ids = data[0].split()
        if mail_ids:
            # Opcional: Marcar los correos como leídos
            for mail_id in mail_ids:
                mail.store(mail_id, '+FLAGS', '\\Seen')
            return True
        return False
    except Exception as e:
        logging.error(f"Error al verificar el correo electrónico para la solicitud: {e}")
        return False

def main():
    # Analizar argumentos de línea de comandos
    parser = argparse.ArgumentParser(description='Script de Monitoreo de Simulaciones de Dinámica Molecular')
    parser.add_argument('--input_dir', type=str, default='.', help='Directorio de entrada que contiene las carpetas de proteínas.')
    parser.add_argument('--output_dir', type=str, default='.', help='Directorio de salida para informes y gráficos.')
    parser.add_argument('--sender_email', type=str, required=True, help='Dirección de correo electrónico del remitente.')
    parser.add_argument('--recipient_email', type=str, required=True, help='Dirección de correo electrónico del destinatario.')
    parser.add_argument('--email_password', type=str, help='Contraseña del correo electrónico del remitente.')
    parser.add_argument('--interval_hours', type=int, default=6, help='Intervalo en horas para el informe periódico.')
    parser.add_argument('--email_subject', type=str, default='Informe MD', help='Línea de asunto para desencadenar la generación del informe.')
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    sender_email = args.sender_email
    recipient_email = args.recipient_email
    email_password = args.email_password
    if email_password is None:
        email_password = getpass.getpass(prompt='Contraseña del correo electrónico: ')
    interval_hours = args.interval_hours
    email_subject = args.email_subject

    # Función para generar el informe y enviar el correo electrónico
    def job():
        generate_monitoring_report_and_plots(input_dir, output_dir)

        # Preparar los detalles del correo electrónico
        subject = 'Informe de Monitoreo de Simulaciones MD'
        body = 'Adjunto encontrará el informe de monitoreo y los gráficos generados.'
        attachments = [
            str(output_dir / 'complete_monitoring_report.txt'),
            str(output_dir / 'global_progress.png'),
            str(output_dir / 'progress_per_protein.png'),
            str(output_dir / 'storage_per_protein.png')
        ]

        send_email_report(sender_email, email_password, recipient_email, subject, body, attachments)
        logging.info(f"Informe de monitoreo y gráficos enviados a {recipient_email}.")

    # Programar la tarea periódica
    schedule.every(interval_hours).hours.do(job)

    # Iniciar el bucle principal
    while True:
        # Verificar si hay solicitud por correo
        if check_email_for_request(sender_email, email_password, email_subject):
            logging.info("Solicitud de informe recibida por correo electrónico.")
            job()
        schedule.run_pending()
        time.sleep(60)  # Esperar 60 segundos antes de verificar nuevamente

if __name__ == "__main__":
    main()
