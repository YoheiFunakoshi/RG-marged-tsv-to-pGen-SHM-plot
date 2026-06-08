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
- `<sample>_pgen_shm_scatter_reads.png`
- `<sample>_pgen_shm_kde_unweighted.png`
- `<sample>_pgen_shm_kde_weighted.png`
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
- `pGen=0` の行は `<sample>_pgen_shm_rows.xlsx` と `<sample>_pgen_shm_points.tsv` には残しますが、scatter/KDE plotからは除外します。

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
