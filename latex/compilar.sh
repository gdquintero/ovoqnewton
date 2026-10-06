#!/bin/bash
# Compila el paper revisado, la version con cambios marcados (latexdiff contra la
# version sometida) y la carta de respuesta. Los PDF quedan en latex/pdf/.
set -e
cd "$(dirname "$0")"
BUILD=$(mktemp -d)
cp acqovo.tex acqovo_sometido.tex respuesta_revisor.tex *.pdf "$BUILD"/ 2>/dev/null || true
cd "$BUILD"

latexdiff acqovo_sometido.tex acqovo.tex > acqovo_diff.tex 2>/dev/null
# Las referencias a Assumption A2 (eliminada) solo aparecen en el texto tachado
sed -i 's/\\begin{document}/\\begin{document}\\expandafter\\def\\csname r@a2\\endcsname{{A2}{}{}{}{}}/' acqovo_diff.tex

for f in acqovo acqovo_diff respuesta_revisor; do
    pdflatex -interaction=nonstopmode $f.tex > /dev/null 2>&1 || true
    pdflatex -interaction=nonstopmode $f.tex > /dev/null 2>&1 || true
    echo "$f: $(grep -c '^!' $f.log) errores, $(grep -c 'undefined' $f.log) referencias indefinidas"
done

mkdir -p "$OLDPWD/pdf"
cp acqovo.pdf "$OLDPWD/pdf/acqovo_revisado.pdf"
cp acqovo_diff.pdf "$OLDPWD/pdf/acqovo_cambios.pdf"
cp respuesta_revisor.pdf "$OLDPWD/pdf/"
grep -E "newlabel\{(a1|a3|a4|lem:step|teo:suff-descent|lem:transfer|teo:complexity)\}" acqovo.aux \
    | sed 's/\\newlabel{\([^}]*\)}{{\([^}]*\)}.*/  \1 -> \2/'
rm -rf "$BUILD"
