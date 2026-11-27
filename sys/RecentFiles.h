#ifndef _RecentFiles_h_
#define _RecentFiles_h_

#include "melder.h"

// Ramyses: implementação de uma lista de arquivos recentes para o Praat, com persistência em disco.

constexpr integer MAX_RECENT_FILES = 10;

// Adiciona um arquivo ao topo da lista de recentes e salva no disco
void RecentFiles_add(MelderFile file);

// Carrega a lista salva ao iniciar o Praat
void RecentFiles_load();

// Salva a lista atual no arquivo de preferências do Praat
void RecentFiles_save();

void RecentFiles_installMenu();
void RecentFiles_refreshMenu();

/*
// Exemplo no final do reconhecimento de arquivo em Data_readFromFile (sys/Data.cpp):
if (object) {
    RecentFiles_add(file);
    return object;
}
*/
#endif