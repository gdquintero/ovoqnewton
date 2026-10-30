#!/bin/bash
# Compila el paper revisado, la version con cambios marcados (latexdiff contra la
# version R1 enviada el 2026-10-12) y la carta de respuesta de la segunda ronda.
# Los PDF quedan en latex/pdf/.
set -e
cd "$(dirname "$0")"
BUILD=$(mktemp -d)
cp acqovo.tex acqovo_R1.tex respuesta_revisor_R2.tex *.pdf "$BUILD"/ 2>/dev/null || true
cd "$BUILD"

latexdiff acqovo_R1.tex acqovo.tex > acqovo_diff.tex 2>/dev/null

for f in acqovo acqovo_diff respuesta_revisor_R2; do
    pdflatex -interaction=nonstopmode $f.tex > /dev/null 2>&1 || true
    pdflatex -interaction=nonstopmode $f.tex > /dev/null 2>&1 || true
    echo "$f: $(grep -c '^!' $f.log) errores, $(grep -c 'undefined' $f.log) referencias indefinidas"
done

mkdir -p "$OLDPWD/pdf"
cp acqovo.pdf "$OLDPWD/pdf/version_final.pdf"
cp acqovo_diff.pdf "$OLDPWD/pdf/version_final_cambios.pdf"
cp respuesta_revisor_R2.pdf "$OLDPWD/pdf/"
grep -E "newlabel\{(cor:adaptive|ccotasub|shift|ovo-armijo|theta-kkt|transfer-const|cota-sigma|cota-fcnt|cota-iter|lem:transfer|teo:complexity|tab:active)\}" acqovo.aux \
    | sed 's/\\newlabel{\([^}]*\)}{{\([^}]*\)}.*/  \1 -> \2/'
rm -rf "$BUILD"
