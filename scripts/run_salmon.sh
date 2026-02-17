#!/bin/bash

# Rutas
INDEX="salmon_index"
FASTQ_DIR="data/raw_fastq"
OUTPUT_DIR="results/salmon_quants"

# Crear carpeta de salida
mkdir -p $OUTPUT_DIR

# Lista de muestras (SRR IDs)
SAMPLES=(
"SRR28830406"
  "SRR28830407"
  "SRR28830408"
  "SRR28830409"
  "SRR28830410"
  "SRR28830411"
)

# Loop para cuantificar cada muestra
for SAMPLE in "${SAMPLES[@]}"
do echo "Processing $SAMPLE..."

salmon quant \
    -i $INDEX \
    -l A \
    -1 ${FASTQ_DIR}/${SAMPLE}_1.fastq \
    -2 ${FASTQ_DIR}/${SAMPLE}_2.fastq \
    -p 4 \
    --validateMappings \
    -o ${OUTPUT_DIR}/${SAMPLE}

echo "$SAMPLE done!"
done

echo "All samples quantified!"


