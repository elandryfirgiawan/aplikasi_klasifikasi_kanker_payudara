from django.shortcuts import render
from .forms import UploadFileForm, UploadFastaForm

import re
from collections import defaultdict, Counter
from Bio import SeqIO
import pandas as pd
import joblib
import matplotlib.pyplot as plt
import os
from sklearn.preprocessing import StandardScaler
from io import BytesIO, TextIOWrapper
import base64

# Load Model
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
model_path = os.path.join(BASE_DIR, 'models')
svm = joblib.load(os.path.join(model_path, 'svm_model.pkl'))
x_mean = joblib.load(os.path.join(model_path, 'x_mean.pkl'))
x_std = joblib.load(os.path.join(model_path, 'x_std.pkl'))
pca = joblib.load(os.path.join(model_path, 'pca_model.pkl'))

# =======================
# FUNGSI PENDUKUNG
# =======================

def build_pattern_stats(sequences):
    pattern_stats = defaultdict(Counter)
    for seq in sequences:
        for i in range(len(seq) - 2):
            left = seq[i]
            middle = seq[i + 1]
            right = seq[i + 2]
            if left in "ATGC" and middle in "ATGC" and right in "ATGC":
                pattern_stats[(left, right)][middle] += 1
    return pattern_stats

def replace_all_N(sequence, pattern_stats):
    sequence = list(sequence)

    # Tahap 1: Tangani 'N' tunggal dahulu
    changed = True
    while changed:
        changed = False
        for i in range(1, len(sequence) - 1):
            if sequence[i] == 'N' and sequence[i - 1] in 'ATGC' and sequence[i + 1] in 'ATGC':
                left = sequence[i - 1]
                right = sequence[i + 1]
                candidates = pattern_stats.get((left, right), {})
                if candidates:
                    replacement = candidates.most_common(1)[0][0]
                    sequence[i] = replacement
                    changed = True

    # Tahap 2: Tangani blok N berurutan
    while 'N' in sequence:
        for i in range(len(sequence)):
            if sequence[i] == 'N':
                start = i
                while i < len(sequence) and sequence[i] == 'N':
                    i += 1
                end = i

                left = sequence[start - 1] if start > 0 else None
                right = sequence[end] if end < len(sequence) else None

                if left in "ATGC" and right in "ATGC":
                    candidates = pattern_stats.get((left, right), {})
                    if candidates:
                        replacement = candidates.most_common(1)[0][0]
                        sequence[start] = replacement
                        break
                break
    return ''.join(sequence)


def build_markov_chain_single_sequence(seq):
    transition_counts = defaultdict(lambda: defaultdict(int))
    total_transitions = defaultdict(int)
    seq = re.sub(r'[^AGCT]', '', seq.upper())
    for i in range(len(seq) - 1):
        current_nuc = seq[i]
        next_nuc = seq[i + 1]
        transition_counts[current_nuc][next_nuc] += 1
        total_transitions[current_nuc] += 1
    feature_vector = {}
    for current_nuc in "AGCT":
        for next_nuc in "AGCT":
            feature_vector[f"{current_nuc}->{next_nuc}"] = transition_counts[current_nuc].get(next_nuc, 0) / total_transitions[current_nuc] if total_transitions[current_nuc] > 0 else 0
    return feature_vector

# =======================
# VIEW INDEX
# =======================

def index(request):
    result = None
    plot_url = None
    fasta_result = None

    file_form = UploadFileForm()
    fasta_form = UploadFastaForm()

    if request.method == 'POST':
        # Excel file
        if 'file' in request.FILES:
            file_form = UploadFileForm(request.POST, request.FILES)
            if file_form.is_valid():
                excel_file = request.FILES['file']
                df = pd.read_excel(excel_file)

                X = df.drop(columns=['file_name'])
                y = df['file_name'].str.strip().replace({'Invasive': 'Invasif', 'Non-Invasive': 'Non-Invasif'})

                Z = (X - x_mean) / x_std
                x_pca = pca.transform(Z)
                pred = svm.predict(x_pca)
                df['Predicted'] = pred
                result = df[['file_name', 'Predicted']].to_html(index=False)

                numerical_labels = y.map({'Invasif': 0, 'Non-Invasif': 1})
                plt.figure(figsize=(6, 4))
                plt.scatter(x_pca[:, 0], x_pca[:, 1], c=numerical_labels, s=20, edgecolors='k')
                plt.title("Visualisasi PCA")
                plt.xlabel("PC1")
                plt.ylabel("PC2")
                buffer = BytesIO()
                plt.savefig(buffer, format='png')
                buffer.seek(0)
                image_png = buffer.getvalue()
                buffer.close()
                plot_url = base64.b64encode(image_png).decode('utf-8')

        # FASTA file (single sequence, header optional)
        elif 'fasta_file' in request.FILES:
            fasta_form = UploadFastaForm(request.POST, request.FILES)
            if fasta_form.is_valid():
                fasta_file = request.FILES['fasta_file']
                text_stream = TextIOWrapper(fasta_file.file, encoding='utf-8')
                lines = text_stream.readlines()

                # Jika ada header
                if lines[0].startswith(">"):
                    try:
                        text_stream.seek(0)
                        record = SeqIO.read(text_stream, "fasta")
                        sequence = str(record.seq).upper()
                    except Exception:
                        fasta_result = "Gagal membaca file FASTA"
                        return render(request, 'breast_cancer_app/index.html', {
                            'file_form': file_form,
                            'fasta_form': fasta_form,
                            'result': result,
                            'fasta_result': fasta_result,
                            'plot_url': plot_url
                        })
                else:
                    # Tidak ada header: gabungkan semua baris
                    sequence = ''.join([line.strip() for line in lines]).upper()

                if 'N' in sequence:
                    pattern_stats = build_pattern_stats([sequence])
                    sequence = replace_all_N(sequence, pattern_stats)

                fitur_dict = build_markov_chain_single_sequence(sequence)
                X_baru = pd.DataFrame([fitur_dict])
                X_baru = X_baru[x_mean.index]
                Z_baru = (X_baru - x_mean) / x_std
                Z_baru = Z_baru.fillna(0)
                X_baru_pca = pca.transform(Z_baru)
                prediksi = svm.predict(X_baru_pca)
                fasta_result = f"{prediksi[0]}"

    return render(request, 'breast_cancer_app/index.html', {
        'file_form': file_form,
        'fasta_form': fasta_form,
        'result': result,
        'fasta_result': fasta_result,
        'plot_url': plot_url
    })
