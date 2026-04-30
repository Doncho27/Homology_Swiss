import os
import requests
import warnings
import gzip
from Bio.PDB import PDBParser, PDBIO
from Bio.Data.IUPACData import protein_letters_3to1
from Bio import Align
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

def extract_sequence_robust(pdb_file):
    """Extrae todos los residuos del PDB (monómero) sin importar el ID de la cadena."""
    parser = PDBParser(QUIET=True)
    
    try:
        with gzip.open(pdb_file, 'rt') as f:
            structure = parser.get_structure('estructura', f)
    except (OSError, EOFError):
        structure = parser.get_structure('estructura', pdb_file)
    
    valid_residues = []
    pdb_seq = ""
    
    # Recorre todas las cadenas en el primer modelo (ideal para monómeros)
    for chain in structure[0]:
        for res in chain.get_residues():
            if res.id[0] == ' ':
                try:
                    letra = protein_letters_3to1[res.resname.capitalize()]
                    pdb_seq += letra
                    valid_residues.append(res)
                except KeyError:
                    pass
                    
    if not pdb_seq:
        raise ValueError("No se encontraron aminoácidos válidos en el archivo.")
        
    return pdb_seq, structure, valid_residues

def get_uniprot_sequence(uniprot_id):
    print(f"[*] Descargando FASTA de UniProt ({uniprot_id}) una sola vez...")
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.strip().split('\n')
        return ''.join(lines[1:])
    else:
        raise ConnectionError("No se pudo descargar la secuencia.")

def renumber_and_save(structure, valid_residues, pdb_seq, uniprot_seq, output_file):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    best_alignment = aligner.align(uniprot_seq, pdb_seq)[0]
    
    mapping = {}
    for (u_start, u_end), (p_start, p_end) in zip(best_alignment.aligned[0], best_alignment.aligned[1]):
        for i in range(u_end - u_start):
            mapping[p_start + i] = u_start + i + 1 

    original_ids = [res.id for res in valid_residues]
    
    for i, res in enumerate(valid_residues):
        res.id = (res.id[0], -10000 - i, res.id[2])
            
    for i, res in enumerate(valid_residues):
        if i in mapping:
            res.id = (res.id[0], mapping[i], res.id[2])
        else:
            try:
                res.id = original_ids[i]
            except ValueError:
                pass

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)
#-------------------------------------- MODIFY HERE --------------------------------------
def main():
    CARPETAS = ["7LHW"]
    UNIPROT_ID = "Q5S007"

    try:
        uniprot_seq = get_uniprot_sequence(UNIPROT_ID)
    except Exception as e:
        print(f"Error fatal: {e}")
        return

    for carpeta in CARPETAS:
        if not os.path.exists(carpeta):
            continue
            
        print(f"\n--- Procesando carpeta: {carpeta} ---")
        for archivo in os.listdir(carpeta):
            if archivo.startswith("modelo_") and archivo.endswith(".pdb") and not archivo.endswith("_rn.pdb"):
                ruta_entrada = os.path.join(carpeta, archivo)
                nombre_base, extension = os.path.splitext(archivo)
                ruta_salida = os.path.join(carpeta, f"{nombre_base}_rn{extension}")
                
                print(f"[*] Reenumerando: {archivo} -> {nombre_base}_rn{extension}")
                
                try:
                    # Ya no le pasamos 'CADENA'
                    pdb_seq, structure_obj, valid_residues = extract_sequence_robust(ruta_entrada)
                    renumber_and_save(structure_obj, valid_residues, pdb_seq, uniprot_seq, ruta_salida)
                except Exception as e:
                    print(f"    Error con {archivo}: {e}")

    print("\nProceso de reenumeracion masiva finalizado.")

if __name__ == "__main__":
    main()
