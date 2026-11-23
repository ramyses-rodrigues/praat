#form Salvar áudio com FFmpeg
#    sentence Pasta_de_destino /Users/seu_usuario/Desktop
#    word Formato mp3
#    comment Exemplos de formato: mp3, ogg, flac, m4a, aac
#    word Bitrate 192k
#endform

if numberOfSelected("Sound") <> 1    
	exitScript: "Por favor, selecione exatamente um objeto Sound na lista."
endif

nome$ = selected$("Sound")
# 2. Abre a janela nativa de "Salvar arquivo"
# O segundo argumento sugere o nome padrão com a extensão desejada (ex: .mp3, .ogg, .m4a)

caminho_final$ = chooseWriteFile$: "Salvar áudio convertido como...", nome$ + ".mp3"

# Interrompe o script silenciosamente se o usuário clicar em "Cancelar"
if caminho_final$ = ""
	exitScript()
endif

nome$ = selected$("Sound")
caminho_wav$ = pasta_de_destino$ + "/" + nome$ + "_temp.wav"
caminho_final$ = pasta_de_destino$ + "/" + nome$ + "." + formato$

# 1. Salva o arquivo temporário em WAV
Save as WAV file: caminho_wav$

# 2. Executa o FFmpeg para converter
# Nota: No Windows, se o ffmpeg não estiver no PATH, coloque o caminho completo (ex: "C:\ffmpeg\bin\ffmpeg.exe")
#comando$ = "ffmpeg -y -i """ + caminho_wav$ + """ -b:a " + bitrate$ + " """ + caminho_final$ + """"
#runSystem: comando$
asynchronous runSubprocess: "ffmpeg", "-y", "-i", caminho_wav$, -b:a, bitrate$, caminho_final$
#appendInfoLine: "Tarefa concluída"

# 3. Remove o arquivo WAV temporário
deleteFile: caminho_wav$

writeInfoLine: "Arquivo salvo com sucesso em: ", caminho_final$
