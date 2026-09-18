# MonSNP-Pipeline：可复现的高密度 SNP 芯片设计流程

[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A57.0-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxxx)

[English README](README.md)

MonSNP-Pipeline 把一组 WGS/PLINK 数据转换成可用于 SNP 芯片设计的候选位点和严格
探针面板。你只需要提供基因型、参考基因组、设计窗口和需要使用的注释或区域文件，
流程会完成频率计算、LD tag set 筛选、位点评分、候选排序和 120-nt 探针筛选。

它有两种常见用途：

- **使用自己的群体或物种数据设计芯片**；
- **复现蒙古马 58K 芯片的论文结果**（需要另外下载论文数据包）。

仓库内置一个 6-SNP 小型数据集，可先用它确认安装和基本流程，再投入真实数据。

## 主要特点

- **从基因型到探针面板**：从 PLINK BED/BIM/FAM 开始，自动完成候选筛选和探针设计。
- **评分方式可配置**：可组合 MAF、VEP、保守区、QTL、表达数据、选择信号或自定义区域。
- **可先小范围试跑**：一条染色体验证通过后，再切换到全基因组。
- **可分阶段检查**：在 LD 计算和探针筛选前后都能停下来检查结果。
- **运行环境可复现**：使用 Snakemake 和 Conda 管理任务及软件版本。

## 你将得到什么

完整运行后，最重要的文件位于 `<results_dir>/06_final_panel/`。其中
`<results_dir>` 是配置文件中的 `results_dir`（主配置默认为 `results`）：

| 文件 | 用途 |
|---|---|
| `final_strict_probe_panel.tsv` | 最终严格筛选的探针集合，包含位点、探针序列、排名和评分 |
| `final_strict_probe_panel.bed` | 最终位点的基因组坐标 |
| `summary_report.txt` | 每个主要筛选阶段的数量汇总 |
| `validation.json` | 序列、坐标和唯一性检查结果 |

如果只运行到评分阶段，候选位点在
`<results_dir>/02_scoring/selection_ld_tag_set_candidates.tsv`，适合先检查排序是否符合预期。

## 5 分钟安装检查

运行环境：Linux x86-64、Conda/Miniforge 和 Git LFS。

```bash
git lfs install
git clone https://github.com/YourName/MonSNP-Pipeline.git
cd MonSNP-Pipeline
conda env create -f environment.yaml
conda activate monsnp-pipeline
bash demo/run_demo.sh
```

看到 `Demo assertions passed` 即表示核心安装检查通过。Demo 会把结果写入
`demo/work/`，并停在 OligoMiner 之前，因此不会触发大型探针比对。Demo 的输入和预期
结果见 [`demo/README.md`](demo/README.md)。

## 用自己的数据运行

### 第一步：确认输入文件

至少需要以下文件：

| 文件 | 说明 |
|---|---|
| `PREFIX.bed/.bim/.fam` | PLINK 二进制基因型文件；三者必须共享同一前缀 |
| `assembly.fa` | 与基因型坐标完全一致的参考基因组 |
| `assembly.fa.fai` | FASTA 索引；没有时可由工具自动生成 |
| `design_windows.bed` | 0-based、half-open 的设计窗口 |

如果使用论文模式，还需要 VEP、保守区、QTL、表达和选择信号等注释文件。只想使用
自己的区域权重时，可以采用 `custom` 模式，不需要论文专用注释轨迹。

### 第二步：准备 PLINK 和设计窗口

根目录下的 `tools/` 是运行主流程前的独立预处理工具，不会原位修改源文件。

```bash
python tools/prepare_plink_input.py \
  --bfile data/raw \
  --reference data/assembly.fa \
  --out data/prepared/population \
  --chr-set 31 \
  --threads 8

python tools/make_windows.py \
  --reference data/assembly.fa \
  --output data/design_windows.bed \
  --window-size 40000 \
  --step 40000
```

马通常使用 `--chr-set 31`；其他物种请改成该物种在 PLINK 中使用的常染色体数。
输入格式和坐标约定见 [`docs/input_formats.md`](docs/input_formats.md)。

### 第三步：复制配置并做预检

自有群体和自定义区域建议从以下模板开始：

```bash
cp config/custom_regions.example.yaml config/my_run.yaml
```

然后至少修改 `inputs.bfile`、`inputs.reference`、`inputs.reference_fai`、
`inputs.windows`，以及 `scoring.weighted_regions` 中的区域路径和权重。完成后运行：

```bash
python tools/preflight.py --config config/my_run.yaml --base-dir . --report preflight.json
```

预检显示 `PASS` 或 `PASS_WITH_WARNINGS` 后再启动 Snakemake。完整字段说明见
[`config/README.md`](config/README.md)；评分规则见
[`docs/custom_scoring.md`](docs/custom_scoring.md)。

该模板默认使用 `results_dir: results_custom`。如果希望结果写入 `results/`，请在配置
中把它改为 `results_dir: results`。

## 推荐的运行方式

真实 WGS 数据不建议第一次就跑全流程。先在一条染色体上试跑，确认输入和中间结果，
再改成全基因组：

```yaml
# 快速试跑
target_chromosomes:
  - ECA1

# 确认无误后再运行全部窗口
# target_chromosomes: all
```

### 三个阶段

| 阶段 | 命令 | 运行到哪里 | 适合什么时候 |
|---|---|---|---|
| A：注释准备 | `run_annotation` | 窗口、频率和可选 VEP 注释 | 第一次检查输入 |
| B：候选排序 | `run_scoring` | LD tag sets、综合评分和候选位点 | 检查排序、进入探针筛选前 |
| C：完整结果 | `run_full` | OligoMiner、严格探针和最终面板 | 确认中间结果后 |

示例命令（将 `config/my_run.yaml` 换成你的配置）：

```bash
# 阶段 A：快速检查输入和注释
snakemake -s workflow/Snakefile \
  --configfile config/my_run.yaml --cores 8 \
  --software-deployment-method conda run_annotation

# 阶段 B：计算 LD 并生成候选位点
snakemake -s workflow/Snakefile \
  --configfile config/my_run.yaml --cores 16 \
  --software-deployment-method conda run_scoring

# 阶段 C：运行完整探针设计
snakemake -s workflow/Snakefile \
  --configfile config/my_run.yaml --cores 32 \
  --software-deployment-method conda run_full
```

> **耗时提醒：** LD tag set 计算和 OligoMiner 探针筛选是流程中最耗时的部分，可能
> 需要大量 CPU、内存和磁盘空间。全基因组运行前请先完成一条染色体的测试，并使用
> `--dry-run` 检查任务数量。`--cores` 应根据服务器或调度系统的实际资源设置。

## 复现论文结果

论文复现需要从 GigaDB/论文数据归档下载完整输入文件，将其放入
`data/manuscript/`，再按 [`docs/manuscript_reproduction.md`](docs/manuscript_reproduction.md)
更新 `config/config.yaml`。之后运行：

```bash
snakemake -s workflow/Snakefile \
  --configfile config/config.yaml --cores 32 \
  --software-deployment-method conda run_full
```

最终面板的数量和 SHA-256 应与 `resources/manuscript_checksums.tsv` 中的记录核对。
正式发布前还需补充真实 DOI、数据登录号、仓库地址和论文引用信息。

## 进一步阅读

- [输入格式和坐标约定](docs/input_formats.md)
- [自有数据快速指南](docs/custom_data_quickstart.zh-CN.md)
- [自定义区域评分](docs/custom_scoring.md)
- [方法和算法](docs/methods.md)
- [版本历史与可复现性](docs/manuscript_reproduction.md)
- [数据可用性](docs/data_availability.md)
- [发布检查清单](docs/release_checklist.md)

## 引用与许可证

如果使用 MonSNP-Pipeline，请引用：

> Wang Y, et al. (2026). *Title of the paper*. GigaScience. DOI: 待补充。

发布前请更新 `CITATION.cff`。本项目采用 [MIT License](LICENSE)，内置 OligoMiner
文件的许可说明见 [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)。
