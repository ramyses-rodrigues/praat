import os
import re

# Configurações do escaneamento
DIRETORIO_PROJETO = "./"  # Varre a partir de onde o script for executado
ARQUIVO_SAIDA = "projeto.dot"

# Regex para capturar tanto #include "modulo.h" quanto #include <modulo.h>
REGEX_INCLUDE = re.compile(r'#include\s*["<]([^">]+)[">]')

# Lista de extensões que queremos analisar
EXTENSOES_CONTEUDO = ('.cpp', '.c', '.cc', '.cxx', '.h', '.hpp', '.hxx')

# Includes comuns da biblioteca padrão (STL) para ignorar (evita poluir o grafo)
STL_IGNORES = {
    'iostream', 'vector', 'string', 'memory', 'map', 'set', 'algorithm', 
    'cmath', 'fstream', 'sstream', 'iomanip', 'functional', 'thread', 'mutex'
}

def extrair_relacionamentos():
    nos_fontes = set()
    nos_headers = set()
    arestas = set()

    for raiz, _, arquivos in os.walk(DIRETORIO_PROJETO):
        for arquivo in arquivos:
            if arquivo.endswith(EXTENSOES_CONTEUDO):
                caminho_completo = os.path.join(raiz, arquivo)
                
                # Classifica o nó para colorir depois
                if arquivo.endswith(('.cpp', '.c', '.cc', '.cxx')):
                    nos_fontes.add(arquivo)
                else:
                    nos_headers.add(arquivo)

                try:
                    with open(caminho_completo, 'r', encoding='utf-8', errors='ignore') as f:
                        for linha in f:
                            match = REGEX_INCLUDE.search(linha)
                            if match:
                                incluido = match.group(1)
                                nome_incluido = os.path.basename(incluido)
                                
                                # Filtra se for biblioteca padrão do C++
                                if nome_incluido not in STL_IGNORES:
                                    arestas.add((arquivo, nome_incluido))
                except Exception as e:
                    print(f"Erro ao ler {arquivo}: {e}")
                    
    return nos_fontes, nos_headers, arestas

def gerar_dot(nos_fontes, nos_headers, arestas):
    with open(ARQUIVO_SAIDA, 'w', encoding='utf-8') as f:
        f.write("digraph G {\n")
        f.write("    rankdir=LR;\n")  # Layout da esquerda para a direita
        f.write("    node [fontname='Helvetica', fontsize=11];\n\n")
        
        # Estiliza os arquivos fontes (.cpp) como retângulos azuis
        f.write("    // Arquivos de Origem\n")
        f.write("    node [shape=box, style=filled, color='#A2C8EC', fontcolor='#1E3A5F'];\n")
        for no in sorted(nos_fontes):
            f.write(f'    "{no}";\n')
            
        # Estiliza os cabeçalhos (.h) como elipses cinzas
        f.write("\n    // Arquivos de Cabeçalho\n")
        f.write("    node [shape=ellipse, style=filled, color='#E0E0E0', fontcolor='#333333'];\n")
        for no in sorted(nos_headers):
            f.write(f'    "{no}";\n')
            
        # Cria as conexões de relacionamento
        f.write("\n    // Relacionamentos (#include)\n")
        f.write("    edge [color='#555555', arrowhead=vee, arrowsize=0.8];\n")
        for origem, destino in sorted(arestas):
            f.write(f'    "{origem}" -> "{destino}";\n')
            
        f.write("}\n")
    print(f" Sucesso! Arquivo '{ARQUIVO_SAIDA}' gerado com sucesso.")

if __name__ == "__main__":
    fontes, headers, conexoes = extrair_relacionamentos()
    gerar_dot(fontes, headers, conexoes)
