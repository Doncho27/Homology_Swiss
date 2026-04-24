import requests
import time
import os
from datetime import datetime
from Bio import SeqIO

# Configuracion
API_TOKEN = "c04c035a5b3d79953dc075bfcb88426a50409a47"
CARPETAS = ["G2019S", "I2020T"]
POST_URL = "https://swissmodel.expasy.org/user_template"
HEADERS = {"Authorization": f"Token {API_TOKEN}"}

def obtener_secuencia(ruta_carpeta):
    for archivo in os.listdir(ruta_carpeta):
        if archivo.endswith((".fasta", ".fa")):
            ruta_completa = os.path.join(ruta_carpeta, archivo)
            registros = list(SeqIO.parse(ruta_completa, "fasta"))
            if registros:
                return str(registros[0].seq)
    return None

def limpiar_pantalla():
    os.system('cls' if os.name == 'nt' else 'clear')

def ejecutar_flujo():
    proyectos = []
    
    print("Iniciando envio por lotes. Esto tomara un momento...")
    
    # Fase 1: Enviar todos
    for mutacion in CARPETAS:
        if not os.path.exists(mutacion):
            continue
            
        secuencia = obtener_secuencia(mutacion)
        if not secuencia:
            continue

        for archivo_pdb in os.listdir(mutacion):
            if archivo_pdb.endswith(".pdb") and not archivo_pdb.startswith("modelo_"):
                ruta_pdb = os.path.join(mutacion, archivo_pdb)
                
                with open(ruta_pdb, 'r') as f_pdb:
                    coordenadas = f_pdb.read()
                
                payload = {
                    "target_sequences": [secuencia],
                    "template_coordinates": coordenadas,
                    "project_title": f"{mutacion}_{archivo_pdb}"
                }
                
                response = requests.post(POST_URL, headers=HEADERS, json=payload)
                
                if response.ok:
                    proyectos.append({
                        'id': response.json()['project_id'],
                        'carpeta': mutacion,
                        'archivo': archivo_pdb,
                        'estado': 'ENVIADO'
                    })
                else:
                    proyectos.append({
                        'id': 'N/A',
                        'carpeta': mutacion,
                        'archivo': archivo_pdb,
                        'estado': 'ERROR_ENVIO'
                    })

    # Fase 2: Monitoreo cada 2 minutos (120 segundos)
    procesos_activos = True
    
    while procesos_activos:
        limpiar_pantalla()
        hora = datetime.now().strftime("%H:%M:%S")
        print(f"--- Tablero de Modelado SWISS-MODEL (Actualizado: {hora}) ---")
        print(f"{'Carpeta':<10} | {'Archivo PDB':<15} | {'Estado Actual':<15}")
        print("-" * 45)
        
        procesos_activos = False
        
        for p in proyectos:
            # Si el proyecto ya termino, no consultamos la API, solo imprimimos su estado
            if p['estado'] not in ['DESCARGADO', 'FAILED', 'ERROR', 'ERROR_ENVIO']:
                procesos_activos = True
                url_estado = f"https://swissmodel.expasy.org/project/{p['id']}/models/summary/"
                
                try:
                    r = requests.get(url_estado, headers=HEADERS)
                    if r.ok:
                        data = r.json()
                        nuevo_estado = data.get('status', p['estado'])
                        
                        if nuevo_estado == 'COMPLETED':
                            modelos = data.get('models', [])
                            if modelos and modelos[0].get('coordinates_url'):
                                url_descarga = modelos[0]['coordinates_url']
                                res_pdb = requests.get(url_descarga, headers=HEADERS)
                                
                                ruta_salida = os.path.join(p['carpeta'], f"modelo_{p['archivo']}")
                                with open(ruta_salida, "wb") as f:
                                    f.write(res_pdb.content)
                                p['estado'] = 'DESCARGADO'
                            else:
                                p['estado'] = 'SIN_COORDENADAS'
                        else:
                            p['estado'] = nuevo_estado
                except Exception:
                    pass # Si hay un error de red puntual, el estado se mantiene igual hasta el proximo intento
            
            # Imprimir la linea de la tabla
            print(f"{p['carpeta']:<10} | {p['archivo']:<15} | {p['estado']:<15}")
        
        if procesos_activos:
            print("\nEsperando 5 minutos para la siguiente actualizacion...")
            time.sleep(120)

    print("\nProceso finalizado. Puedes revisar las carpetas.")

if __name__ == "__main__":
    ejecutar_flujo()
