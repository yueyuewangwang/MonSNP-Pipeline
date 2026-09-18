# MonSNP-Pipeline 目录结构

本文档说明公开仓库中的源码目录和运行时输出。所有路径均相对于仓库根目录。

## 源码目录树

```text
MonSNP-Pipeline/
├── .github/workflows/
├── config/
├── demo/
│   └── input/
├── docs/
├── resources/
├── tests/
│   ├── data/
│   │   └── plink/
│   ├── expected/
│   └── mock_bin/
├── tools/
├── workflow/
│   ├── envs/
│   ├── schemas/
│   ├── scripts/
│   └── vendor/OligoMiner/
├── .gitattributes
├── .gitignore
├── .snakemake-workflow-catalog.yml
├── CHANGELOG.md
├── CITATION.cff
├── CONTRIBUTING.md
├── DIRECTORY_STRUCTURE.md
├── LICENSE
├── Makefile
├── README.md
├── README.zh-CN.md
├── THIRD_PARTY_NOTICES.md
└── environment.yaml
```

## 主要目录

### `config/`

- `config.yaml`：蒙古马论文分析配置模板，包含完整注释轨迹和正式流程参数。
- `custom_regions.example.yaml`：自有群体与自定义加权区域配置示例。
- `manuscript_bundle_manifest.example.tsv`：论文归档数据包清单示例。
- `README.md`：配置键、染色体范围、LD 参数和阶段目标说明。

配置中的 `target_chromosomes` 可设为 `all` 或显式染色体列表。实际运行时建议先复制
模板为 `config/my_run.yaml`。

### `demo/`

包含 6 个 SNP、4 个样本、1-kb FASTA 和两个加权区域的小型演示。
`run_demo.sh` 依次执行 PLINK 转换、输入规范化、窗口生成、预检、LD tag sets 和评分，
`assert_demo.py` 核对关键结果。生成文件写入 `demo/work/`，不属于源码。

### `docs/`

包括自有数据快速指南、输入格式、自定义评分、算法说明、版本历史与可复现性、数据
可用性、软件引用和发布检查清单。

### `resources/`

- `scoring_weights.tsv`：可配置且可追踪的评分权重。
- `manuscript_checksums.tsv`：论文冻结输入与结果的 SHA-256 记录。

### `tests/`

- `data/`：合成 FASTA、FAI、BED、VEP、AFREQ、表达和注释轨迹输入。
- `data/plink/`：真实格式的小型 PLINK BED/BIM/FAM 测试集。
- `mock_bin/`：仅用于无外部依赖测试的轻量 PLINK/PLINK2 模拟程序。
- `config*.yaml`：预计算输入及 PLINK 输入两条测试路径。
- `assert_*.py` 与 `run_smoke.sh`：数值、格式和关键逻辑测试。

### `tools/`

该目录保留在仓库根目录，并明确定义为 **Pre-processing utilities（独立预处理工具）**。
它们在 Snakemake 主流程前运行：

- `prepare_plink_input.py`：规范 PLINK REF/ALT 和 `CHROM:POS:REF:ALT` ID。
- `make_windows.py`：生成 FASTA FAI 和 BED4 设计窗口。
- `preflight.py`：检查文件、坐标、参考等位基因、染色体范围、权重及依赖。

这些工具不会原位修改源文件，均支持 `--help`。

### `workflow/`

- `Snakefile`：主工作流和 `run_annotation`、`run_scoring`、`run_full` 三个阶段目标。
- `envs/`：核心 Python、PLINK 和 OligoMiner 规则环境。
- `schemas/config.schema.yaml`：启动时执行的配置结构校验。
- `scripts/chromosome_partitions.py`：创建窗口、AFREQ 和 LD tag-set 染色体分片。
- `scripts/build_oligominer_resources.sh`：一次性建立共享 Bowtie2/Jellyfish 资源。
- `scripts/` 其余文件：VEP 评分、候选排序、探针筛选、面板组装和结果验证。
- `vendor/OligoMiner/`：冻结的 OligoMiner 源码与原始许可证。

LD 和 OligoMiner 任务在 `by_chrom/` 中分片执行，随后由 gather 规则生成规范汇总路径。

## 根目录文件

- `README.md`、`README.zh-CN.md`：中英文项目入口。
- `environment.yaml`：Demo 和工作流管理环境。
- `Makefile`：测试、预检、dry-run 和三个阶段目标的快捷命令。
- `.gitattributes`：为 BED/BIM/FAM/VCF.GZ/BCF/PGEN/BAM/CRAM/JF 等格式启用 Git LFS。
- `.gitignore`：排除结果、缓存、本地配置和大型临时文件。
- `CITATION.cff`、`LICENSE`、`THIRD_PARTY_NOTICES.md`：引用与许可信息。

## 运行时输出

```text
results/
├── 00_inputs/
│   ├── by_chrom/        # 目标窗口分片
│   └── target_windows.bed
├── 00_plink/
│   ├── by_chrom/        # AFREQ 与 LD tag-set 分片
│   ├── variants.afreq
│   └── window_ld_tag_sets.tsv
├── 01_annotation/       # VEP 功能注释评分
├── 02_scoring/          # 分项总分和每窗口候选排序
├── 03_probe_candidates/
│   └── by_chrom/        # 参考校验后的候选侧翼分片
├── 04_oligominer/
│   ├── shared/          # 共享 Bowtie2/Jellyfish 资源
│   └── by_chrom/        # 染色体筛选结果
├── 05_probe_selection/  # 汇总后的严格探针选择
├── 06_final_panel/      # 最终面板、汇总报告和验证结果
└── logs/

.snakemake/              # Snakemake 状态与规则环境缓存
demo/work/               # Demo 生成文件
tests/work/              # 测试生成文件
```

研究级输出和大型输入不应直接提交 Git；应存入稳定数据仓库并保存 SHA-256 校验和。
