# Gera as dependências entre módulos Fortran para o Makefile.
#
# Uso: awk -v objdir=DIR -f mkdeps.awk pass=1 ARQUIVOS... pass=2 ARQUIVOS...
#
# Passo 1: registra em qual arquivo cada MODULE é definido.
# Passo 2: para cada USE, emite "objdir/arquivo.o: objdir/modulo_definido_em.o".
# Módulos externos (OMP_LIB, ISO_*, ...) são ignorados porque não estão em
# nenhum dos arquivos listados.

function objname(path,    n, parts, base) {
    n = split(path, parts, "/")
    base = parts[n]
    sub(/\.[fF]90$/, ".o", base)
    return objdir "/" base
}

{
    line = tolower($0)
    sub(/\r$/, "", line)          # fontes com terminação CRLF
    sub(/!.*/, "", line)
}

pass == 1 && line ~ /^[ \t]*module[ \t]+[a-z0-9_]+[ \t]*$/ {
    name = line
    sub(/^[ \t]*module[ \t]+/, "", name)
    sub(/[ \t]*$/, "", name)
    if (name != "procedure") defined[name] = FILENAME
    next
}

pass == 2 && line ~ /^[ \t]*use[ \t,]/ {
    name = line
    sub(/^[ \t]*use[ \t]*(,[ \t]*(non_)?intrinsic[ \t]*::)?[ \t]*/, "", name)
    sub(/[^a-z0-9_].*$/, "", name)
    if ((name in defined) && defined[name] != FILENAME) {
        key = FILENAME SUBSEP name
        if (!(key in seen)) {
            seen[key] = 1
            print objname(FILENAME) ": " objname(defined[name])
        }
    }
}
