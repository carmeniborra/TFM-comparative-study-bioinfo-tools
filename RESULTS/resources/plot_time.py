import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def extract_metrics(filenames):
    """
    Extrae walltime y memoria de los timelines de Nextflow.
    """
    task_regex = re.compile(r'\{"label": "(?P<label>[^"]*)",.*?"times": \[(?P<times>.*?)\]\}')
    label_regex = re.compile(r'"label": "(?P<dur_str>[\d\.dhms\s]+) \\/ (?P<mem_str>[\d\. \w]+)"')

    data = []
    for fname in filenames:
        try:
            with open(fname, 'r', encoding='utf-8') as f:
                content = f.read()
            
            for match in task_regex.finditer(content):
                process_name = match.group('label')
                times_block = match.group('times')
                
                dur_matches = label_regex.finditer(times_block)
                for d_match in dur_matches:
                    d_str = d_match.group('dur_str')
                    m_str = d_match.group('mem_str')
                    
                    # 1. Convertir Tiempo a Minutos
                    total_ms = 0
                    d = re.search(r'(\d+)d', d_str)
                    h = re.search(r'(\d+)h', d_str)
                    m = re.search(r'(\d+)m', d_str)
                    s = re.search(r'(\d+)s', d_str)
                    if d: total_ms += int(d.group(1)) * 86400000
                    if h: total_ms += int(h.group(1)) * 3600000
                    if m: total_ms += int(m.group(1)) * 60000
                    if s: total_ms += int(s.group(1)) * 1000
                    
                    # 2. Convertir Memoria a GB
                    mem_gb = 0
                    mem_val = re.search(r'(\d+\.?\d*)', m_str)
                    if mem_val:
                        val = float(mem_val.group(1))
                        if 'GB' in m_str: mem_gb = val
                        elif 'MB' in m_str: mem_gb = val / 1024
                    
                    # 3. Mapear Herramienta
                    tool = None
                    if "KRAKEN2_KRAKEN2" in process_name: tool = "Kraken2"
                    elif "BRACKEN" in process_name: tool = "Bracken"
                    elif "MOTUS_PROFILE" in process_name: tool = "mOTUs"
                    elif "METAPHLAN_METAPHLAN" in process_name: tool = "MetaPhlAn"
                    elif "KAIJU_KAIJU" in process_name: tool = "Kaiju"
                    elif "GANON_CLASSIFY" in process_name: tool = "Ganon"
                    elif "CENTRIFUGE_CENTRIFUGE" in process_name: tool = "Centrifuge"
                    elif "DIAMOND_BLASTX" in process_name: tool = "DIAMOND"
                    
                    if tool and total_ms > 0:
                        data.append({'Herramienta': tool, 'Minutos': total_ms / 60000.0, 'RAM_GB': mem_gb})
        except Exception as e:
            print(f"Error en {fname}: {e}")
            
    return pd.DataFrame(data)

# --- LISTA COMPLETA DE TUS ARCHIVOS ---
archivos = [
    'execution_timeline_2025-11-20_11-53-04.html', 
    'execution_timeline_2025-11-27_16-40-58.html',
    'execution_timeline_2025-11-30_15-15-14.html',
    'execution_timeline_2025-12-04_10-11-22.html',
    'execution_timeline_2025-12-05_13-38-47.html', 
    'execution_timeline_2025-12-09_10-46-07.html',
    'execution_timeline_2025-12-13_21-04-49.html', 
    'execution_timeline_2025-12-17_09-09-22.html'
]

df = extract_metrics(archivos)
orden = ['Bracken', 'Kraken2', 'mOTUs', 'MetaPhlAn', 'Kaiju', 'Ganon', 'Centrifuge', 'DIAMOND']
df['Herramienta'] = pd.Categorical(df['Herramienta'], categories=orden, ordered=True)

# --- GRÁFICO DE TIEMPOS ---
plt.figure(figsize=(12, 6))
sns.set_style("whitegrid")
ax = sns.boxplot(x='Herramienta', y='Minutos', data=df, palette='Spectral')
ax.set_yscale('log')
plt.yticks([1, 10, 60, 600, 1440, 5000], ['1m', '10m', '1h', '10h', '1d', '3.5d'])
plt.title('Walltime por Muestra (Consolidado)', fontsize=14)
plt.ylabel('Tiempo (Escala Log)')
plt.savefig('boxplot_tiempos_final.png', dpi=300)
plt.show()

# --- GRÁFICO DE MEMORIA (OPCIONAL PARA 3.2.2) ---
plt.figure(figsize=(12, 6))
ax_mem = sns.boxplot(x='Herramienta', y='RAM_GB', data=df, palette='viridis')
plt.title('Consumo de Memoria RAM (Max RSS)', fontsize=14)
plt.ylabel('Memoria (GB)')
plt.savefig('boxplot_memoria_final.png', dpi=300)
plt.show()