
# Abre janela do Windows para selecionar o arquivo
arquivo$ = chooseReadFile$: "Selecione o arquivo de áudio ou vídeo"

# Encerra o script caso o usuário clique em Cancelar
if arquivo$ = ""    
	exitScript: "Operação cancelada ou arquivo não selecionado."
endif

# Define o caminho do arquivo WAV temporário na pasta temporária do sistema
arquivo_wav$ = temporaryDirectory$ + "/praat_ffmpeg_temp.wav"

# Converte para WAV (mono, 44100 Hz, PCM 16-bit) usando o FFmpeg
# As aspas triplas (""") protegem caminhos com espaços no Windows
writeInfoLine: "Iníciando tarefa. Aguarde..."
asynchronous runSubprocess: "ffmpeg", "-y", "-i", arquivo$, "-vn", "-ac", "1", "-ar", "44100", "-c:a", "pcm_s16le", arquivo_wav$
appendInfoLine: "Tarefa concluída"
# runSystem: "ffmpeg -y -i """, arquivo$, """ -vn -ac 1 -ar 44100 -c:a pcm_s16le """, arquivo_wav$, """"

# Carrega o WAV convertido na janela de objetos do Praat
som = Read from file: arquivo_wav$

# Extrai apenas o nome original do arquivo (sem pasta e extensão) para renomear o objeto
nome_base$ = replace_regex$ (arquivo$, "^.*[\\/]", "", 0)
nome_base$ = replace_regex$ (nome_base$, "\.[^.]+$", "", 0)
Rename: nome_base$

# Remove o arquivo WAV temporário do disco após carregar na memória
deleteFile: arquivo_wav$