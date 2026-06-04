#!/bin/bash
set -e

# ==============================================================================
# --- НАСТРОЙКИ ----------------------------------------------------------------
# ==============================================================================
INPUT_DIR="./astral_pipeline/perfect_aligned_genes"
PARALLEL_JOBS=3       # Сколько генов считаем одновременно
THREADS_PER_JOB=6     # Сколько ядер на один ген
# ==============================================================================

# Считаем точное количество исходных файлов
TOTAL_FILES=$(ls "$INPUT_DIR"/*.fasta 2>/dev/null | wc -l)

if [ "$TOTAL_FILES" -eq 0 ]; then
    echo "Ошибка: В папке $INPUT_DIR нет файлов .fasta!"
    exit 1
fi

echo "=== Запуск массового построения генных деревьев (IQ-TREE) ==="
echo "Всего генов: $TOTAL_FILES | Ядер задействовано: $((PARALLEL_JOBS * THREADS_PER_JOB))"
echo "Расчет запущен. Пожалуйста, подождите..."
echo "-------------------------------------------------------------"

# 1. Функция для одного гена
run_iqtree() {
    local aln_file="$1"
    local threads="$2"
iqtree -s "$aln_file" -st DNA -m GTR+G -B 1000 -T "$threads" --quiet -redo
}
export -f run_iqtree


# 2. ФОНОВЫЙ СЧЕТЧИК (будет работать параллельно с расчетами)
(
    while true; do
        sleep 10 # Обновляем каждые 10 секунд
        done_count=$(ls "$INPUT_DIR"/*.treefile 2>/dev/null | wc -l)
        awk -v d="$done_count" -v t="$TOTAL_FILES" 'BEGIN { printf "⏳ Построено деревьев: %d из %d | Прогресс: %.2f%%\n", d, t, (d/t)*100 }'
    done
) & 
# Запоминаем ID (PID) нашего фонового счетчика
MONITOR_PID=$!


# 3. ОСНОВНОЙ ЗАПУСК IQ-TREE
# parallel работает без --bar, чтобы не мусорить на экране
find "$INPUT_DIR" -maxdepth 1 -name "*.fasta" | parallel -j "$PARALLEL_JOBS" run_iqtree {} "$THREADS_PER_JOB"
# 4. ЗАВЕРШЕНИЕ
# Как только parallel закончил работу, останавливаем фоновый счетчик
kill $MONITOR_PID 2>/dev/null
wait $MONITOR_PID 2>/dev/null || true

echo "-------------------------------------------------------------"
echo "✅ ПРОЦЕСС ЗАВЕРШЕН!"
echo "Все $TOTAL_FILES генных деревьев успешно построены."



