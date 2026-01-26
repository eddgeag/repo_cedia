#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
echo "Uso: ./run.sh DX043-25 DX044-25 DX045-25 ..."
exit 1
fi

RSCRIPT_BIN=$(which Rscript)
MAIN_SCRIPT="main.R"
SLEEP_TIME=10
LOG_DIR="./logs"

mkdir -p "$LOG_DIR"

if [ ! -x "$RSCRIPT_BIN" ]; then
echo "ERROR: Rscript no encontrado en PATH"
exit 1
fi

if [ ! -f "$MAIN_SCRIPT" ]; then
echo "ERROR: No existe $MAIN_SCRIPT"
exit 1
fi

echo ">>> Iniciando corrida secuencial de exomas"
echo ">>> Logs en: $LOG_DIR"
echo ">>> Muestras: $*"
echo

for muestra in "$@"; do
LOG_FILE="${LOG_DIR}/${muestra}.log"

echo "===================================================="
echo ">>> Procesando muestra: $muestra"
echo ">>> Log: $LOG_FILE"
echo ">>> Inicio: $(date)"

{
  echo ">>> ===== START ${muestra} ====="
  echo ">>> Fecha inicio: $(date)"
  echo
  
  "$RSCRIPT_BIN" "$MAIN_SCRIPT" "$muestra"
  
  echo
  echo ">>> Fecha fin: $(date)"
  echo ">>> ===== END ${muestra} ====="
} &> "$LOG_FILE"

echo ">>> Finalizado: $muestra"
echo ">>> Fin: $(date)"
echo ">>> Sleep ${SLEEP_TIME}s antes de la siguiente muestra"
sleep "$SLEEP_TIME"
done

echo "===================================================="
echo ">>> Todas las muestras finalizadas correctamente"
