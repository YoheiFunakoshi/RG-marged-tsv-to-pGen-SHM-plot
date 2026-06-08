# RG marged AIRR TSV to pGen-SHM plot

IgBLAST `-outfmt 19` の AIRR TSV から、pGen、SHM、pGen-SHM plot を作成するWindows向けGUIです。

このソフトの入力は FASTA ではありません。前段の merged FASTA -> IgBLAST AIRR TSV 変換を済ませた `.igblast.airr.tsv` を入力します。

## 出力

入力TSVを1つ選んで実行すると、出力フォルダに以下を作成します。

- `<sample>_qc_summary.tsv`
- `<sample>_pgen_bins.tsv`
- `<sample>_pgen_bins_unique.png`
- `<sample>_pgen_bins_weighted.png`
- `<sample>_shm_hist.tsv`
- `<sample>_shm_hist.png`
- `<sample>_shm_hist_weighted.tsv`
- `<sample>_shm_hist_weighted.png`
- `<sample>_pgen_shm_rows.xlsx`
- `<sample>_pgen_shm_points.tsv`
- `<sample>_pgen_shm_beta1_unique_junction_points.tsv`
- `<sample>_pgen_shm_roi_summary.tsv`
- `<sample>_pgen_shm_scatter_reads.png`
- `<sample>_pgen_shm_kde_unweighted.png`
- `<sample>_pgen_shm_kde_log_density.png`
- `<sample>_pgen_shm_kde_weighted.png`
- `<sample>_pgen_shm_kde_weighted_log_density.png`
- `<sample>_pgen_shm_kde_beta1_unique_junction_unweighted.png`
- `<sample>_pgen_shm_kde_beta1_unique_junction_weighted.png`
- `<sample>_pgen_shm_scatter_weighted.png`
- `<sample>_run_log.txt`
- `pgen_cache.tsv`

## Word仕様書

処理単位ごとの詳細仕様は、以下のWord文書に分けています。

- `TSVからSHMデータ作成仕様書_20260607.docx`
- `TSVからpGenデータ作成仕様書_20260607.docx`
- `TSVからpGen-SHMデータ作成仕様書_20260607.docx`

## beta_1 解析ルール

- 対象は基本的に IGH です。
- `locus` 列があり、値が `IGH` 以外の行は除外します。
- `locus` が空欄の行は保持し、QC summaryに件数を記録します。
- `productive == T` を採用します。
- `vj_in_frame` 列がある場合は `T` を採用します。
- `stop_codon` 列がある場合は `F` を採用します。
- `junction` は A/C/G/T のみを採用します。
- `junction_aa` は空欄、`*`、`X` を除外します。
- SHMは `v_identity` から計算します。
  - 0-100スケール: `100 - v_identity`
  - 0-1スケール: `(1 - v_identity) * 100`
- pGenは OLGA `human_B_heavy` の `compute_aa_CDR3_pgen(junction_aa)` で計算します。
- 採用されたAIRR各行を中間データとして `<sample>_pgen_shm_rows.xlsx` に出力します。
- 中間Excelは採用行のみを含み、フィルタで除外された行は含めません。
- 中間Excelの基本列は `sequence_id`, `junction_aa`, `shm`, `pgen`, `log10_pgen`, `junction`, `v_identity`, `locus`, `productive`, `v_call`, `j_call`, `same_xy_count`, `plot_weight` です。
- pGen-SHM plotの1点は、基本的にフィルタ後のAIRR 1行、すなわち1つの `sequence_id` です。
- pGen-SHM plotのx軸はその行の `junction_aa` の `log10(pGen)`、y軸はその行の `v_identity` から計算したSHMです。
- 同一の `log10_pgen` と `shm` を持つ行が複数ある場合、Excelでは各 `sequence_id` 行を保持し、描画時だけ `same_xy_count` / `plot_weight` を点サイズおよびKDE重みに反映します。
- `weighted` 出力は、通常は同一 `junction` / `junction_aa` がTSV内に出た行数で重み付けします。
- `duplicate_count` 列があり、GUIで `Use duplicate_count for weighted outputs if present` にチェックした場合は、weighted出力の重みに `duplicate_count` を使います。
- SHM histogramとpGen-SHM KDEは、非重み付き版と重み付き版の両方を出力します。
- `scatter_reads` は点の存在確認、通常KDEは高密度な主集団の確認、`kde_log_density` は低密度だが意味のある集団の確認に使います。
- `beta1_unique_junction` 出力は旧Beta1解析との比較用です。1点を unique `junction` とし、SHMは同一junction内のmedianを使います。
- `pgen_shm_roi_summary.tsv` は低SHM・高pGen領域などの固定ROIについて、row-level方式とBeta1互換方式の件数を並べて出力します。
- `pGen=0` の行は `<sample>_pgen_shm_rows.xlsx` と `<sample>_pgen_shm_points.tsv` には残しますが、scatter/KDE plotからは除外します。

## pGen-SHM解析の考え方

pGen-SHM plotは探索的な新しい表示であり、現時点では唯一の正解図を決めません。
このGUIでは、元データを潰さない row-level Excel/TSV を中心データとして保存し、そこから複数の見方を出します。

- 主解析単位は、フィルタ後のAIRR 1行、すなわち原則として1つの `sequence_id` です。
- `junction_aa` からpGenを計算し、同じ行の `v_identity` からSHMを計算します。
- 同じ `junction_aa` でもSHMが異なる可能性があるため、主解析ではmedianに潰さず各点を残します。
- Beta1互換出力では、旧解析との比較用に unique `junction` ごとのSHM medianを1点として保存します。
- 図の見え方だけで判断せず、scatter、log-density KDE、ROI summaryを併用します。

## pGen-SHM図の使い分け

- `*_pgen_shm_scatter_reads.png`: 点の存在確認用です。低密度集団が本当に存在するかを見るのに向きます。
- `*_pgen_shm_kde_unweighted.png`: row-levelのlinear KDEです。read密度の主集団を確認します。
- `*_pgen_shm_kde_weighted.png`: row-level weighted KDEです。観測read量を反映するため、生物学的な主図候補です。
- `*_pgen_shm_kde_log_density.png`: 低密度だが意味のある集団を見落とさないための確認図です。
- `*_pgen_shm_kde_beta1_unique_junction_unweighted.png`: 前任者法・Beta1との比較用です。多様性寄りの見え方になります。
- `*_pgen_shm_roi_summary.tsv`: 低SHM・高pGen領域などを、図の色ではなく数値で確認するための表です。

## 議論点

- row-level方式は情報を潰さない一方で、read数、PCR増幅、samplingの影響を受けます。
- Beta1互換方式はunique junction diversityを見やすい一方で、同一junction内のSHMばらつきをmedianに潰します。
- linear KDEでは高密度集団が強く出るため、低密度集団が見えにくくなることがあります。
- log-density KDEは低密度集団の確認に有用ですが、主図としては広がりすぎる場合があります。
- 生物学的な存在量を押す場合は row-level weighted KDE を主図候補にし、scatterとROI summaryで補助確認する方針が妥当です。

## 使い方

1. `run_gui.bat` をダブルクリックします。
2. `AIRR TSV` に IgBLAST outfmt 19 のTSVを指定します。
3. `Output folder` を確認します。
4. `Run pGen + SHM + pGen-SHM plot` を押します。

最初に `Check setup` を押すと、`numpy`, `matplotlib`, `scipy`, `olga` と OLGA `human_B_heavy` モデルの有無を確認できます。

`duplicate_count` 列を重みとして使いたい場合だけ、`Use duplicate_count for weighted outputs if present` にチェックします。通常のAIRR TSVに `duplicate_count` 列がない場合は、TSV内の行数が重みになります。

`Recalculate all pGen (ignore existing cache)` はデフォルトONです。この場合、既存の `pgen_cache.tsv` は読み込まず、毎回OLGAで全unique `junction_aa` のpGenを再計算します。計算が正常終了した後、新しい `pgen_cache.tsv` に置き換えます。

`pGen workers` はOLGA pGen計算の並列数です。通常作業をしながら使う場合は、動作確認済みの `6` を推奨します。PCを解析専用にできる場合だけ、必要に応じて増やします。

## 依存関係

このフォルダ直下の `.venv` に入れる想定です。

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install numpy pandas matplotlib scipy olga openpyxl
```

## コマンドライン実行

GUIを使わずにCLIからも実行できます。

```powershell
.\.venv\Scripts\python.exe .\airr_pgen_shm_plot_beta1.py --input sample.igblast.airr.tsv --outdir result --sample SAMPLE --pgen-workers 6
```

CLIもデフォルトではpGenを全再計算します。既存cacheを使いたい場合だけ `--use-pgen-cache` を付けます。

## 注意

- 生成物や検体データはGitHubへ置かない運用です。
- `pgen_cache.tsv` は `junction_aa -> pGen` のキャッシュです。デフォルト設定では既存cacheを信用せず、毎回再計算してからcacheを作り直します。
- V領域が短いreadはSHM推定が不安定になり得ます。必要に応じてGUIの `Min V alignment length` を設定します。
