# dl_binder_design 설치 순서 (Colab)

RFdiffusion 백본에 서열을 올리고 계면을 평가하는 다음 단계입니다.
`github.com/nrbennet/dl_binder_design` — RFdiffusion 논문과 같은 그룹의 표준 후속 파이프라인.

설치 전에 읽을 것: 아래 ①~③이 모두 통과해야 ④로 갑니다.
②에서 막히면 "가벼운 경로"(맨 아래)로 우회하십시오.

---

## 무엇을 하는 도구인가

| 단계 | 하는 일 | 자원 |
|---|---|---|
| `silentfrompdbs` | PDB 묶음 → silent 파일 하나 | CPU, 즉시 |
| `proteinmpnn_fastrelax` | binder 서열 설계 + FastRelax | **CPU**, PyRosetta 필요 |
| `af2_initial_guess` | 복합체 재예측, `pae_interaction` 산출 | **GPU**, JAX |

`pae_interaction < 10`이 보고된 in-silico 성공 기준이고, 사전 등록 §5의 H1이
이 값을 씁니다.

---

## ⓿ 게이트가 두 개입니다 — 어느 쪽을 어디서 도는가

H1 은 세 단계이고, 지표가 쓸 만한지 검정하는 게이트가 두 군데에 있습니다.
이 문서는 Colab 의 1·2 단계만 다룹니다. 3 단계와 그 게이트는 로컬입니다.

| | 단계 | 어디서 | 게이트 |
|---|---|---|---|
| 1 | ProteinMPNN (서열 설계) | Colab, CPU | — |
| 2 | AF2 initial guess | Colab, **GPU** | **B.** native PD-1 의 `pae_interaction` < 10 |
| 3 | PyRosetta 계면 점수 | **로컬 WSL2** | **A.** 4ZQK 복합체의 `dG_separated` < 0 |

A 를 먼저 도십시오. 1 분이고, 여기서 막히면 Colab 에 1·2 를 설치할
이유가 줄어듭니다. (2026-11-04 통과: `dG_separated` -32.74 REU)

1·2 단계는 노트북으로 정리해 두었습니다.
`notebooks/H1_dl_binder_design.ipynb`. 아래 ①~③ 은 그 노트북이 무엇을
왜 하는지에 대한 배경이고, 실행은 노트북에서 하십시오. 노트북은
`.ipynb` 를 손으로 고치지 않기 위해 `build_h1_notebook.py` 로 생성합니다.

```bash
cd /mnt/c/Users/PC/Downloads/pdl1-binder-design &&   /mnt/c/Users/PC/Downloads/exercise/.venv-linux/bin/python -I   scripts/h1_rosetta.py --positive-control
```

Ala 스캔과 같은 venv 를 직접 호출합니다(레포 1 의 `.venv-linux`).
GPU 불필요, Colab 과 무관합니다. 스크립트는 `dG_separated` 가 음수가
아니면 exit 1 로 멈춥니다. native 복합체가 통과하지 못하는 기준으로
설계를 판정할 수 없습니다.

B 는 2 단계를 돌릴 때 설계군보다 먼저 native PD-1(4ZQK 사슬 B)로
도십시오. 이전 작업에서 AF2 의 ipTM 이 native 에서 0.12 가 나와 계면에
대해 아무 말도 할 수 없게 된 적이 있고, 그때는 양성 대조군을 늦게 돌려서
늦게 알았습니다.

둘 다 실패하면 H1 은 측정 불가이고, 그것을 결론으로 보고합니다.
한쪽만 통과하면 통과한 지표로만 H1 을 보고하고 나머지를 한계에 적습니다.

---

## ⓵ 실측 설치 기록 (2026-11-04) — 15개 실패 지점

아래 ②~④의 계획과 충돌하면 이 절을 따르십시오. 계획은 실행 전에 쓴
것이고 이 절은 실제로 통과한 경로입니다. 설치에만 반나절이 걸렸고 실패가
15번 났습니다. 그중 4번은 제(어시스턴트) 판단 착오였습니다. 표에
그대로 남깁니다.

재현용 노트북은 `notebooks/H1_dl_binder_design.ipynb`이고, 손으로 고치지
말고 생성기 `notebooks/build_h1_notebook.py`를 고쳐 다시 만드십시오
(.ipynb는 JSON이라 손편집하면 깨집니다). 아래 수정은 전부 생성기에
반영되어 있습니다.

| # | 증상 | 원인 | 조치 |
|---|---|---|---|
| 1 | `silentfrompdbs: not found` | `silent_tools`가 레포 최상위가 아니라 `include/silent_tools`에 있음 | 경로를 하드코딩하지 말고 **자동 탐색** + `silentfrompdbs` 존재까지 확인 |
| 2 | `conda: command not found` | 노트북을 **다시 업로드하면 새 세션**이고 `/content`가 비워짐 | 중간 셀부터 재개 불가. conda 존재를 검사하는 가드 셀을 둠 |
| 3 | `PackagesNotFoundError: pyrosetta` | 레포의 yml에 RosettaCommons 채널이 없음 | ↓ 4번과 함께 |
| 4 | `conda: error: unrecognized arguments: -c` | conda 26.3.2의 `env create`에 `-c` 없음 | yml을 복사해 `channels` 머리에 `https://conda.graylab.jhu.edu`를 **`yaml.safe_dump`로 주입** (pip 하위 리스트가 보존됨) |
| 5 | `numpy has no attribute 'trapz'` | jax 0.4.20 ↔ numpy 2 충돌 | **내가 `numpy<2`로 내렸는데 방향이 틀렸습니다.** 6번을 유발 |
| 6 | `jax.core has no attribute 'take_current_trace'` | haiku가 jax 0.4.20에 비해 너무 새로움. `ml-dtypes 0.3.2`는 tensorflow 2.21과 충돌 | ↓ 7번에서 일괄 해결 |
| 7 | `devices [CpuDevice(id=0)]`, "CUDA-enabled jaxlib is not installed" | yml의 `jax[cuda12_pip]==0.4.20` 핀이 유지되지 않음 | **내려가지 말고 올라갑니다** — `jax[cuda12]==0.4.34 jaxlib==0.4.34 ml-dtypes>=0.5.1 numpy>=2 dm-haiku==0.0.11` → `[CudaDevice(id=0)]` |
| 8 | `cannot import name 'SCOPData' from 'Bio.Data'` | biopython에서 제거됨. **내가 "1.84에서 제거"라고 했는데 1.83에도 없어** 핀이 실패 | 벤더링된 `alphafold/data/mmcif_parsing.py`에 try/except shim: `SCOPData.protein_letters_3to1` → `PDBData.protein_letters_3to1_extended`. **소스 수정이므로 기록이 남아야 합니다** |
| 9 | ProteinMPNN 가중치 0개 발견 | 이 버전에서 서브모듈이 아님(`.gitmodules` 없음) | `mpnn_fr/` 안에서 `git clone https://github.com/dauparas/ProteinMPNN.git` (`.pt` 12개, `v_48_020.pt` 확인) |
| 10 | `silent_tools`에서 `No module named 'numpy'` | 셰뱅이 `#!/usr/bin/env python`이라 conda base를 집음 | `PATH`에 `/usr/local/envs/<env>/bin`을 **앞에** 붙임 |
| 11 | `cannot import name 'silent_tools' from 'silent_tools'` | **내가 `PYTHONPATH=ST`로 줘서** 패키지가 동명 모듈에 가려짐 | `PYTHONPATH=os.path.dirname(ST)` |
| 12 | `Failed to initialize NumPy: _ARRAY_API not found` | MPNN 환경의 torch 2.0이 numpy 1.x ABI로 빌드됨 | **MPNN 환경에만** `numpy<2` (AF2와 반대) → numpy 1.26.4 + torch 2.0.1 |
| 13 | `has already been processed. Skipping` | 남아 있는 `check.point` | 삭제하고 실행별 디렉터리에 `-checkpoint_name` 지정. RFdiffusion의 cautious 모드와 **같은 함정** |
| 14 | 합친 silent이 0 바이트 | **`silentcat`은 존재하지 않습니다 — 제가 지어냈습니다** | 헤더 두 줄(`SEQUENCE:`, `description`으로 끝나는 `SCORE:`)을 2번째 조각부터 떼고 이어 붙이는 함수를 직접 씀 |
| 15 | 태그 수가 56/57, 63/64처럼 어긋남 | **내가 stdout+stderr를 섞어 읽어** `silentls`의 DeprecationWarning 단어가 태그로 세어짐 | `r.stdout`만 읽음 |

### 환경별 numpy가 서로 반대입니다

이게 가장 함정입니다. 두 환경을 같은 핀으로 맞추려 하면 둘 중 하나가 깨집니다.

| 환경 | numpy | 핵심 핀 |
|---|---|---|
| `af2_binder_design` | **2.4.6** (`numpy>=2`) | `jax[cuda12]==0.4.34`, `jaxlib==0.4.34`, `ml-dtypes 0.6.0`, `dm-haiku==0.0.11` |
| `proteinmpnn_binder_design` | **1.26.4** (`numpy<2`) | `torch 2.0.1` |

### 실측 시간 (Colab L4 23 GB, 12 vCPU, 52 GB RAM)

| 단계 | 시간 |
|---|---|
| conda 환경 2개 생성 | 15분 43초 |
| ProteinMPNN 직렬 | 1.35분/설계 → 50건 **134분** |
| ProteinMPNN 병렬 (12조각, 6워커, `CUDA_VISIBLE_DEVICES=''`) | **25.7분** (5배) |
| AF2 initial guess | 50건 **59분** (1.2분/건) |
| FastRelax (MPNN 내부) | ~70초/설계, **단일 스레드** |

MPNN의 FastRelax가 단일 스레드이므로 병렬화 이득이 거기서 나옵니다. GPU는
오히려 비워야 합니다(`CUDA_VISIBLE_DEVICES=''`) — 6개 프로세스가 한 GPU를
다투면 느려집니다.

### 병렬 단계 뒤에는 반드시 태그를 세십시오

MPNN이 입력 50건 중 3건을 조용히 떨어뜨렸습니다(`pdl1_epi2_14`,
`pdl1_epi_10`, `pdl1_ctrl2_10`). 오류를 내지 않고 출력 48·49건이 되었습니다.
`-debug`로 재실행하니 3/3 전부 "reported success" — 일시적 실패이고 원인은
특정하지 못했습니다.

세지 않았으면 48 대 49로 AF2를 2시간 돌리고 사전 등록한 n이 깨진 것을
몰랐을 것입니다.

```python
def tags_of(silent):
    # stdout 만 읽습니다. stderr 를 섞으면 경고 단어가 태그로 세어집니다 (실패 15번).
    r = subprocess.run(f'{ST}/silentls {silent}', shell=True,
                       capture_output=True, text=True, env=BASE)
    return [t for t in r.stdout.split() if t]
```

### 게이트 B의 결론 — 이 지표는 native 단백질에 못 씁니다

설치를 끝내고 돌린 양성 대조군이 실패했습니다. `pae_interaction` 26.06,
사슬 순서를 바꿔도 27.67. `plddt_target`은 92인데 `plddt_binder`만 33입니다.

원인은 AF2 initial guess가 MSA를 쓰지 않는다는 것입니다
(`predict.py`: `max_extra_msa = 5`, `max_msa_clusters = 5`,
`make_msa_features(msas=[[feat_holder.seq]])`). native Ig β-sandwich는 단일
서열로 접히지 않습니다. 반면 de novo 설계는 이상화된 나선 다발이라 단일
서열로도 접힙니다. 그래서 설계에는 쓸 수 있고 native에는 못 쓰는 지표이고,
결과적으로 양성 대조군으로 native를 고른 것이 잘못된 선택이었습니다.

사전 등록 §4가 이 분기를 미리 규정해 둔 덕에 판정은 막히지 않았습니다.
`dG_separated`로 승격했고 그쪽 양성 대조군은 통과했습니다(−32.74 REU).

---

## ① 두 가지 제약을 먼저 확인

PyRosetta 라이선스. `proteinmpnn_fastrelax`는 PyRosetta를 씁니다.
레포 1에서 내린 판단(공개 교육 자료 = 비상업적)이 그대로 적용되며,
`work/05_rosetta/README_rosetta_note.md`의 고지를 레포 2에도 복사해야 합니다.
이 판단은 본 저자의 것이며 독자에게 이전되지 않습니다.

Colab과 conda. 이 레포는 conda 환경 두 개를 전제합니다. Colab에는 conda가
없으므로 `condacolab`이나 micromamba를 먼저 깔아야 하고, `condacolab.install()`은
런타임을 재시작합니다. 재시작하면 RFdiffusion 쪽 환경(`/content/rfenv`)과
`LD_LIBRARY_PATH` 수정이 모두 사라집니다.

> 따라서 RFdiffusion 배치를 전부 끝내고 Drive에 백업한 뒤
> 별도 세션에서 이 단계를 시작하십시오. 한 세션에서 둘 다 하려 하면
> 런타임 재시작으로 앞 작업을 잃습니다.

---

## ② 환경 구축

```bash
pip -q install condacolab
python -c "import condacolab; condacolab.install()"   # ← 런타임 재시작됨
```

재시작 후:

```bash
git clone https://github.com/nrbennet/dl_binder_design.git /content/dlbd
cd /content/dlbd && git submodule update --init --recursive   # silent_tools
```

환경 두 개를 만듭니다. 레포의 `include/` 안에 환경 yml이 있습니다
(파일명은 버전에 따라 다르므로 `ls /content/dlbd/include/` 로 먼저 확인).

```bash
ls /content/dlbd/include/
mamba env create -f /content/dlbd/include/proteinmpnn_fastrelax.yml
mamba env create -f /content/dlbd/include/af2_binder_design.yml
```

PyRosetta는 PyPI에 없고 전용 채널에서 받습니다.

```bash
mamba install -n proteinmpnn_binder_design \
  -c https://conda.graylab.jhu.edu pyrosetta
```

채널 접근이 막히면 `pyrosetta-installer`(레포 1에서 쓴 방식)로 대체할 수
있는지 확인하십시오 — 같은 바이너리를 가져옵니다.

## ③ 가중치 두 벌

```bash
# ProteinMPNN (서브모듈에 포함되지 않았으면)
ls /content/dlbd/mpnn_fr/ProteinMPNN/vanilla_proteinmpnn/vanilla_model_weights/

# AF2 — params_model_1_ptm.npz 만 필요
mkdir -p /content/dlbd/af2_initial_guess/model_weights/params
cd /content/dlbd/af2_initial_guess/model_weights/params
wget -q https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
tar -xf alphafold_params_2022-12-06.tar params_model_1_ptm.npz
rm alphafold_params_2022-12-06.tar
ls -lh
```

`params_model_1_ptm.npz`가 보이면 통과입니다.

## ④ 실행

```bash
cd /content/dlbd

# PDB -> silent
./include/silent_tools/silentfrompdbs /content/designs/pdl1_*.pdb > all.silent

# 서열 설계 (CPU). 먼저 2건으로 시간을 재십시오
conda run -n proteinmpnn_binder_design python mpnn_fr/dl_interface_design.py \
  -silent all.silent -outsilent mpnn.silent \
  -relax_cycles 1 -seqs_per_struct 1

# 계면 평가 (GPU)
conda run -n af2_binder_design python af2_initial_guess/predict.py \
  -silent mpnn.silent -outsilent af2.silent
```

점수는 `out.sc` (또는 `-scorefilename`으로 지정)에 나오고
`pae_interaction` 열이 H1의 판정값입니다.

비용을 먼저 재십시오. Colab 무료 티어는 vCPU 2개입니다.
FastRelax는 설계당 수십 초~수 분이며, 50건 × 1서열이면 1~3시간일 수 있습니다.
2건으로 시간을 재고 총량을 계산한 뒤 본실행하십시오 — RFdiffusion에서
54분/설계를 한참 들여다본 것과 같은 실수를 반복하지 않기 위해서입니다.

---

## 가벼운 경로 (②에서 막히면)

conda·PyRosetta가 Colab에서 해결되지 않으면, H1을 포기하지 않고도
아래로 갈 수 있습니다. 대신 FastRelax가 빠지므로 보고 시 그 차이를
명시해야 합니다(원 파이프라인과 숫자를 직접 비교하면 안 됩니다).

1. ProteinMPNN 단독 — `github.com/dauparas/ProteinMPNN`, pip만으로 동작.
   `--pdb_path_chains A`(binder)와 `--fixed_positions`로 타깃 고정.
   레포 1에서 이미 쓴 도구이므로 설정을 그대로 재사용할 수 있습니다.
2. AF2 initial guess 대신 ColabDesign AF2-multimer — `ipTM`/`pAE`를 얻습니다.
   단 레포 1에서 native ipTM 0.12가 나왔던 그 경로이므로,
   사전 등록 §4의 양성 대조군(native PD-1)을 반드시 먼저 통과시켜야
   합니다. 통과하지 못하면 H1은 "측정 실패"로 보고하고 H2만 씁니다.
3. PyRosetta 계면 점수 — WSL2 로컬에서 돌립니다(레포 1의 델타 A 스테이지를
   그대로 이식). GPU가 필요 없고 이미 검증된 환경입니다.

---

## H2는 지금 당장 계산할 수 있습니다

사전 등록 §5의 H2는 이 설치를 전혀 기다리지 않습니다. RFdiffusion 출력
PDB만으로 계산됩니다.

```bash
python scripts/h2_hotspot_contacts.py "outputs/pdl1_epi_*.pdb" \
    --hotspots A56,A115,A121,A123 --group epitope
python scripts/h2_hotspot_contacts.py "outputs/pdl1_ctrl_*.pdb" \
    --hotspots A32,A106,A109,A132 --group control
```

H2가 주 지표이므로, 설치가 길어져도 핵심 결과는 먼저 손에 들어옵니다.
