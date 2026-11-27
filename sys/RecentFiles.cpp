#include "RecentFiles.h"
#include "praatP.h"
#include "praat_script.h" // Declara praat_executeScriptFromFile
#include <vector>
#include <string>

static std::vector<std::u32string> s_recentPaths;
static std::vector<std::u32string> s_activeMenuTitles;
static bool s_loaded = false;
static bool s_submenuCreated = false;

static void getRecentStorageFile(MelderFile outFile) {
    MelderFolder prefDir = Melder_preferencesFolder();
    MelderFolder_getFile(prefDir, U"recent_files.txt", outFile);
}

// Escapa aspas duplas (") para ("") conforme a sintaxe de script do Praat
static std::u32string escapePraatQuotes(const std::u32string& text) {
    std::u32string escaped;
    escaped.reserve(text.size() + 8);
    for (char32 ch : text) {
        if (ch == U'"') {
            escaped += U"\"\"";
        } else {
            escaped.push_back(ch);
        }
    }
    return escaped;
}

// Extrai o caminho do arquivo de dentro de uma linha AddedMenuCommand
static std::u32string extractPathFromAddedCommand(const std::u32string& line) {
    const std::u32string marker = U"Read from file: \"\"";
    size_t start = line.find(marker);
    if (start == std::u32string::npos) return U"";

    start += marker.size();
    size_t end = line.rfind(U"\"\"\"");
    if (end == std::u32string::npos || end <= start) return U"";

    std::u32string raw = line.substr(start, end - start);
    std::u32string unescaped;
    for (size_t i = 0; i < raw.size(); i++) {
        if (raw[i] == U'"' && i + 1 < raw.size() && raw[i + 1] == U'"') {
            unescaped.push_back(U'"');
            i++;
        } else {
            unescaped.push_back(raw[i]);
        }
    }
    return unescaped;
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
                    std::u32string path = extractPathFromAddedCommand(line);
                    if (!path.empty()) {
                        s_recentPaths.push_back(path);
                        if ((integer)s_recentPaths.size() >= MAX_RECENT_FILES) break;
                    }
                    line.clear();
                }
            } else {
                line.push_back(*p);
            }
            p++;
        }
        if (!line.empty() && (integer)s_recentPaths.size() < MAX_RECENT_FILES) {
            std::u32string path = extractPathFromAddedCommand(line);
            if (!path.empty()) {
                s_recentPaths.push_back(path);
            }
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

        for (integer i = 0; i < (integer)s_recentPaths.size(); i++) {
            const std::u32string& fullPath = s_recentPaths[i];

            structMelderFile tempFile{};
            Melder_pathToFile(fullPath.c_str(), &tempFile);
            std::u32string fileName = MelderFile_name(&tempFile);

            std::u32string escapedTitle = escapePraatQuotes(fileName);
            std::u32string escapedPath  = escapePraatQuotes(fullPath);

            // Grava no padrão de script AddedMenuCommand do Praat
            MelderString_append(
                &outText,
                U"AddedMenuCommand: \"Objects\", \"Open\", \"",
                i + 1, U". ", escapedTitle.c_str(),
                U"\", \"Open recent\", 1, \"Read from file: \"\"",
                escapedPath.c_str(),
                U"\"\"\"\n"
            );
        }

        MelderFile_writeText_e(&file, outText.string, kMelder_textOutputEncoding::UTF8);
    } catch (MelderError) {
        Melder_clearError();
    }
}

// Mantém a lista dos títulos atualmente exibidos no submenu "Open recent"
static void clearRecentSubmenu() {
    for (const auto& title : s_activeMenuTitles) {
        try {
            // Escapa eventuais aspas no título do arquivo
            std::u32string escapedTitle = escapePraatQuotes(title);

            // Executa o comando nativo de script do Praat para remover o item do menu:
            // Remove menu command: "Objects", "Open", "1. arquivo.wav"
            autoMelderString removeScript;
            MelderString_copy(
                &removeScript,
                U"Remove menu command: \"Objects\", \"Open\", \"",
                escapedTitle.c_str(),
                U"\""
            );

            praat_executeScriptFromText(removeScript.string);
        } catch (MelderError) {
            // Ignora caso o item já tenha sido removido manualmente
            Melder_clearError();
        }
    }
    s_activeMenuTitles.clear();
}


void RecentFiles_refreshMenu() {
    if (s_recentPaths.empty()) return;

    // Limpa todos os itens antigos do submenu antes de atualizar
    clearRecentSubmenu();

    // 1. Cria o submenu pai "Open recent" apenas na primeira execução
    if (!s_submenuCreated) {
        praat_addMenuCommand(
            U"Objects",
            U"Open",
            U"Open recent",
            U"Read from file: ",
            0,
            nullptr
        );
        s_submenuCreated = true;
    }

    // 2. Adiciona cada item diretamente via API C++ (sem chamar praat_executeScriptFromFileNameWithArguments)
    for (integer i = 0; i < (integer)s_recentPaths.size(); i++) {
        const std::u32string& fullPath = s_recentPaths[i];

        structMelderFile tempFile{};
        Melder_pathToFile(fullPath.c_str(), &tempFile);

        autoMelderString title;
        MelderString_copy(&title, i + 1, U". ", MelderFile_name(&tempFile));

        // Evita adicionar o mesmo rótulo duas vezes na mesma sessão
        bool existsInMenu = false;
        for (const auto& existing : s_activeMenuTitles) {
            if (existing == title.string) {
                existsInMenu = true;
                break;
            }
        }
        if (existsInMenu) continue;

        autoMelderString scriptCmd;
        MelderString_copy(&scriptCmd, U"Read from file: \"", fullPath.c_str(), U"\"");

        praat_addMenuCommandScript(
            U"Objects",
            U"Open",
            title.string,
            U"Abrir recent",
            1,
            scriptCmd.string
        );

        s_activeMenuTitles.push_back(title.string);
    }
}


void RecentFiles_refreshMenu2() {
    // 1. Remove os itens anteriores da sessão para evitar duplicação no submenu
    for (const auto& oldTitle : s_activeMenuTitles) {
        try {
            autoMelderString scriptRemove;
            // Monta o comando de script nativo do Praat para remover o item do menu
            MelderString_copy(
                &scriptRemove,
                U"Remove menu command: \"Objects\", \"Open\", \"",
                oldTitle.c_str(),
                U"\""
            );
            praat_executeScriptFromText(scriptRemove.string);
        } catch (MelderError) {
            Melder_clearError();
        }
    }
    s_activeMenuTitles.clear();

    if (s_recentPaths.empty()) return;

    // 2. Cria o submenu pai "Open recent" apenas na primeira execução
    if (!s_submenuCreated) {
        try {
            praat_addMenuCommand(
                U"Objects",
                U"Open",
                U"Open recent",
                U"Read from file: ",
                0,
                nullptr
            );
            s_submenuCreated = true;
        } catch (MelderError) {
            Melder_clearError();
        }
    }

    // 3. Registra os títulos atuais para controle de remoção na próxima atualização
    for (integer i = 0; i < (integer)s_recentPaths.size(); i++) {
        structMelderFile tempFile{};
        Melder_pathToFile(s_recentPaths[i].c_str(), &tempFile);

        autoMelderString title;
        MelderString_copy(&title, i + 1, U". ", MelderFile_name(&tempFile));
        s_activeMenuTitles.push_back(title.string);
    }

    // 4. Executa o recent_files.txt como script do Praat para popular os itens no menu
    structMelderFile file{};
    getRecentStorageFile(&file);
    if (MelderFile_exists(&file)) {
        try {
            autoMelderString scriptCall;
            // Envolve o caminho completo do arquivo em aspas duplas caso haja espaços (ex.: "C:\Users\...")
            MelderString_copy(&scriptCall, U"\"", file.path, U"\"");
            praat_executeScriptFromFileNameWithArguments(scriptCall.string);
            
        } catch (MelderError) {
            Melder_clearError();
        }
    }
}

void RecentFiles_add(MelderFile file) {
    if (!file || MelderFile_isNull(file)) return;

    if (!s_loaded) {
        RecentFiles_load();
    }

    conststring32 fullPath = file->path;
    if (!fullPath || fullPath[0] == U'\0') return;

    std::u32string newPath(fullPath);

    // Move o arquivo para o topo se já existir na lista
    for (auto it = s_recentPaths.begin(); it != s_recentPaths.end(); ) {
        if (*it == newPath) {
            it = s_recentPaths.erase(it);
        } else {
            ++it;
        }
    }

    s_recentPaths.insert(s_recentPaths.begin(), newPath);

    if ((integer)s_recentPaths.size() > MAX_RECENT_FILES) {
        s_recentPaths.resize(MAX_RECENT_FILES);
    }

    RecentFiles_save();
    RecentFiles_refreshMenu();
}

void RecentFiles_installMenu() {
    RecentFiles_load();
    RecentFiles_refreshMenu();
}
