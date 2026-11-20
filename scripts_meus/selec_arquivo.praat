form Abrir arquivo de áudio/vídeo
    folder Pasta C:\Users\ramys\Desktop
    choice Filtro 1
        button Todos os suportados (mp4, mkv, m4a, ogg, opus, wav)
        button Apenas Vídeos (mp4, mkv, mov, avi, webm)
        button Apenas Áudios (m4a, ogg, opus, aac, wma, mp3)
endform

# Lista todos os arquivos da pasta escolhida
lista = Create Strings as file list: "arquivos", pasta$ + "/*"
total = Get number of strings

# Monta a expressão regular de acordo com a escolha
if filtro = 1
    regex$ = "\.(mp4|mkv|mov|avi|webm|m4a|ogg|opus|aac|wma|mp3|wav)$"
elif filtro = 2
    regex$ = "\.(mp4|mkv|mov|avi|webm)$"
else
    regex$ = "\.(m4a|ogg|opus|aac|wma|mp3|wav)$"
endif

# Filtra a lista mantendo apenas as extensões escolhidas
for i from total to 1
    selectObject: lista
    nome$ = Get string: i
    if index_regex (nome$, "(?i)" + regex$) = 0
        Remove string: i
    endif
endfor

# Exibe uma janela suspensa (beginPause) apenas com os arquivos filtrados
selectObject: lista
qtd_filtrada = Get number of strings
if qtd_filtrada = 0
    removeObject: lista
    exitScript: "Nenhum arquivo compatível encontrado nessa pasta."
endif

beginPause: "Selecione o arquivo"
    optionMenu: "Arquivo_escolhido", 1
    for i to qtd_filtrada
        nome$ = Get string: i
        option: nome$
    endfor
clicked = endPause: "Cancelar", "Abrir", 2, 1

if clicked = 1
    removeObject: lista
    exitScript ()
endif

arquivo$ = pasta$ + "/" + arquivo_escolhido$
removeObject: lista

arquivo$ = chooseReadFile$: "Selecione o arquivo de áudio ou vídeo"
if arquivo$ = ""
    exitScript ()
endif

# Verifica se a extensão está na lista permitida (ignora maiúsculas/minúsculas com (?i))
if index_regex (arquivo$, "(?i)\.(mp4|mkv|mov|avi|webm|m4a|ogg|opus|aac|wma|mp3|wav|flac)$") = 0
    exitScript: "Formato não suportado: ", arquivo$
endif

exitScript ()