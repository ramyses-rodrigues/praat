// Ramyses: implementação de uma lista de arquivos recentes para o Praat, com persistência em disco.

#include "RecentFiles.h"
#include <vector>
#include <string>

static std::vector<std::u32string> s_recentPaths;
static bool s_loaded = false;

// Obtém o caminho do arquivo de histórico na pasta de preferências do Praat
// (No Windows: C:\Users\<usuario>\Praat\recent_files.txt)
// C:\Users\<usuario>\AppData\Roaming\Praat
static void getRecentStorageFile(MelderFile outFile) {
    MelderFolder prefDir = Melder_preferencesFolder();
    MelderFolder_getFile(prefDir, U"recent_files.txt", outFile);
    
    // structMelderFolder structprefDir{};
    // MelderFolder prefDir = &structprefDir;
    // prefDir = Melder_preferencesFolder();
    // // Melder_preferencesFolder(&prefDir);
    // MelderFolder_getFile(prefDir, U"recent_files.txt", outFile);
}

void RecentFiles_load() {
    s_recentPaths.clear();
    s_loaded = true;

    structMelderFile file{};
    getRecentStorageFile(&file);

    if (!MelderFile_exists(&file)) {
        return;
    }

    try {
        autostring32 content = MelderFile_readText(&file);
        const char32* p = content.get();
        if (!p) return;

        std::u32string line;
        while (*p != U'\0') {
            if (*p == U'\r' || *p == U'\n') {
                if (!line.empty()) {
                    s_recentPaths.push_back(line);
                    line.clear();
                    if ((integer)s_recentPaths.size() >= MAX_RECENT_FILES) break;
                }
            } else {
                line.push_back(*p);
            }
            p++;
        }
        if (!line.empty() && (integer)s_recentPaths.size() < MAX_RECENT_FILES) {
            s_recentPaths.push_back(line);
        }
    } catch (MelderError) {
        Melder_clearError();
    }
}

void RecentFiles_save() {
    try {
        structMelderFile file{};
        getRecentStorageFile(&file);

        autoMelderString outText;
        for (const auto& path : s_recentPaths) {
            MelderString_append(&outText, path.c_str(), U"\n");
        }

        MelderFile_writeText_e(&file, outText.string, kMelder_textOutputEncoding::UTF8);
    } catch (MelderError) {
        Melder_clearError();
    }
}

void RecentFiles_add(MelderFile file) {
    if (!file || MelderFile_isNull(file)) return;

    if (!s_loaded) {
        RecentFiles_load();
    }

    conststring32 fullPath = file ? file->path : U"";;
    if (!fullPath || fullPath[0] == U'\0') return;

    std::u32string newPath(fullPath);

    // Remove se já existir na lista (para recolocá-lo no topo sem duplicar)
    for (auto it = s_recentPaths.begin(); it != s_recentPaths.end(); ) {
        if (*it == newPath) {
            it = s_recentPaths.erase(it);
        } else {
            ++it;
        }
    }

    // Insere no início (posição 0 = mais recente)
    s_recentPaths.insert(s_recentPaths.begin(), newPath);

    // Mantém apenas os N arquivos mais recentes
    if ((integer)s_recentPaths.size() > MAX_RECENT_FILES) {
        s_recentPaths.resize(MAX_RECENT_FILES);
    }

    RecentFiles_save();
}

integer RecentFiles_count() {
    if (!s_loaded) RecentFiles_load();
    return static_cast<integer>(s_recentPaths.size());
}

conststring32 RecentFiles_getPath(integer index) {
    if (!s_loaded) RecentFiles_load();
    if (index < 0 || index >= (integer)s_recentPaths.size()) {
        return U"";
    }
    return s_recentPaths[index].c_str();
}
