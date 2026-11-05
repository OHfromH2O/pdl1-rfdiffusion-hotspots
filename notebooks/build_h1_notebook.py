#!/usr/bin/env python3
"""H1_dl_binder_design.ipynb 를 생성한다.

왜 생성 스크립트인가 — .ipynb 는 JSON 이고 손으로 고치면 깨집니다.
셀을 파이썬 리스트로 두면 diff 가 읽히고, 수정이 안전합니다.

    python notebooks/build_h1_notebook.py
"""
import io
import json
import os

MD, CODE = "markdown", "code"
CELLS = []


def md(text):
    CELLS.append((MD, text))


def code(text):
    CELLS.append((CODE, text))


# ──────────────────────────────────────────────────────────────────────
md("""# H1 — 계면 품질 (dl_binder_design)

RFdiffusion 백본 100건에 **서열을 올리고**(ProteinMPNN) **복합체를 다시
예측해**(AF2 initial guess) 계면의 질을 잽니다. 사전 등록 §5 의 H1 입니다.

| 단계 | 하는 일 | 자원 |
|---|---|---|
| `silentfrompdbs` | PDB 묶음 → silent 파일 | CPU, 즉시 |
| ProteinMPNN | binder 서열 설계 (타깃 고정) | CPU |
| AF2 initial guess | 복합체 재예측 → `pae_interaction` | **GPU** |

**`pae_interaction < 10`** 이 보고된 in-silico 성공 기준입니다.

---

## 읽고 시작하십시오

**① 이 노트북은 RFdiffusion 세션과 같이 쓸 수 없습니다.**
아래 `condacolab.install()` 이 **런타임을 재시작**하므로 `/content/rfenv` 와
`LD_LIBRARY_PATH` 수정이 사라집니다. RFdiffusion 배치를 모두 끝내고 Drive 에
백업한 뒤 **새 세션**에서 여십시오. (생성은 2026-11-03 에 100/100 완료)

**② 게이트가 먼저입니다.** 설계군 100 건을 돌리기 전에 native PD-1 로
`pae_interaction` 이 작동하는지 확인합니다(셀 8). 이전 작업에서 AF2 의 ipTM 이
native 복합체에서 0.12 가 나와 계면에 대해 아무 말도 할 수 없게 된 적이 있고,
그때는 양성 대조군을 늦게 돌려서 늦게 알았습니다.

**③ 사슬 순서가 중요합니다.** `dl_binder_design` 은 **사슬 A = binder,
사슬 B = 타깃**을 가정합니다. RFdiffusion 출력은 이미 그 순서입니다.
4ZQK 복합체는 반대이므로 `inputs/PD1_PDL1_binderfirst.pdb` 로 미리
뒤집어 두었습니다 — 게이트에는 그 파일을 씁니다.

**④ 각 셀은 다음 셀로 가기 전에 스스로 검사합니다.** `[ok]` 가 없으면
멈추고 출력을 확인하십시오. 조용히 틀린 상태로 진행하는 것이 가장 비쌉니다.
""")

# ── 1 ────────────────────────────────────────────────────────────────
md("""## 1. 환경 사실 수집

무엇이 주어졌는지 먼저 적습니다. GPU 종류와 메모리는 AF2 가 돌지 말지를
결정하고, 나중에 "왜 느렸는가"를 되짚을 때 유일한 근거가 됩니다.""")

code("""import subprocess, sys, platform

def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return (r.stdout + r.stderr).strip()

print("python :", sys.version.split()[0], "/", platform.platform())
print("nvidia-smi:")
print(sh("nvidia-smi --query-gpu=name,memory.total,driver_version "
         "--format=csv,noheader") or "  (GPU 없음)")
print()
free = sh("df -h /content | tail -1")
print("디스크   :", free)
print()
print("[확인] GPU 가 없으면 Colab 메뉴 > 런타임 > 런타임 유형 변경 에서")
print("       T4 이상을 고르고 이 셀을 다시 돌리십시오. AF2 는 GPU 가 필요합니다.")""")

# ── 2 ────────────────────────────────────────────────────────────────
md("""## 2. Drive 연결 + 입력 확인

RFdiffusion 출력 100 건과 게이트용 양성 대조군이 Drive 에 있어야 합니다.
없으면 여기서 멈춥니다 — 뒤 단계를 설치할 이유가 없습니다.""")

code("""import os, glob
from google.colab import drive
drive.mount('/content/drive')

ROOT = '/content/drive/MyDrive/RFdiffusion'      # 필요하면 바꾸십시오
OUT  = f'{ROOT}/outputs'
INP  = f'{ROOT}/inputs'

epi  = sorted(glob.glob(f'{OUT}/pdl1_epi_*.pdb')) + \\
       sorted(glob.glob(f'{OUT}/pdl1_epi2_*.pdb'))
ctrl = sorted(glob.glob(f'{OUT}/pdl1_ctrl_*.pdb')) + \\
       sorted(glob.glob(f'{OUT}/pdl1_ctrl2_*.pdb'))
gate = f'{INP}/PD1_PDL1_binderfirst.pdb'

print(f'설계군 (epitope) : {len(epi)}건')
print(f'대조군 (control) : {len(ctrl)}건')
print(f'게이트용 양성 대조군: {os.path.isfile(gate)}  {gate}')
print()

ok = True
if len(epi) != 50 or len(ctrl) != 50:
    print('[중단] 각 군 50건이어야 합니다. Drive 경로와 파일명을 확인하십시오.')
    ok = False
if not os.path.isfile(gate):
    print('[중단] PD1_PDL1_binderfirst.pdb 가 없습니다.')
    print('       로컬 inputs/ 에서 Drive 의 inputs/ 로 올리십시오.')
    print('       (사슬 A = PD-1 binder, B = PD-L1 타깃으로 뒤집어 둔 파일)')
    ok = False

# 사슬 순서를 실제로 확인합니다. A = binder 가 아니면 뒤 단계가 조용히 틀립니다.
def chains(p):
    c = {}
    for ln in open(p):
        if ln.startswith('ATOM'):
            c[ln[21]] = c.get(ln[21], 0) + 1
    return c

if epi:
    print('설계 예시 사슬:', chains(epi[0]), '<- A가 binder(약 70-100잔기)여야 함')
if os.path.isfile(gate):
    print('게이트 사슬    :', chains(gate), "<- A=776(PD-1), B=874(PD-L1)")
print()
print('[ok] 다음 셀로' if ok else '[중단] 위 메시지를 해결한 뒤 다시 돌리십시오')""")

# ── 3 ────────────────────────────────────────────────────────────────
md("""## 3. conda 설치 — **이 셀은 런타임을 재시작합니다**

`dl_binder_design` 은 conda 환경 두 개를 전제합니다. Colab 에는 conda 가
없으므로 `condacolab` 을 깝니다.

> **재시작 후에는 셀 1·2 를 다시 돌리지 않아도 됩니다**(Drive 는 셀 4 가
> 다시 마운트합니다). 셀 4 로 바로 가십시오.
>
> **주의 — 노트북을 다시 업로드하면 새 세션입니다.** `/content` 가 비워지고
> conda 도 사라지므로 **셀 3 부터** 다시 시작해야 합니다. `condacolab` 의
> 재시작은 커널만 재시작하므로 `/content` 가 보존되지만, 업로드는 다릅니다.
> 실제로 한 번 혼동해 `conda: command not found` 를 만났습니다.""")

code("""!pip -q install condacolab
import condacolab
condacolab.install()        # ← 여기서 런타임이 재시작됩니다. 정상입니다.""")

# ── 4 ────────────────────────────────────────────────────────────────
md("""## 4. (재시작 후) 저장소 받기 + **구조 확인**

환경 yml 의 파일명은 버전마다 다릅니다. 추측하지 않고 **먼저 출력을
봅니다.** 다음 셀에서 그 이름을 적어 넣습니다.""")

code("""import os, subprocess
from google.colab import drive
if not os.path.ismount('/content/drive'):
    drive.mount('/content/drive')

def sh(c):
    print(f'$ {c}')
    r = subprocess.run(c, shell=True, capture_output=True, text=True)
    print((r.stdout + r.stderr).strip()[:4000])
    return r.returncode

# conda 가 없으면 이 세션에서 셀 3 을 아직 돌리지 않은 것입니다.
# 노트북을 **다시 업로드하면 새 세션**이 되어 /content 가 비워지고 conda 도
# 사라집니다. condacolab.install() 의 재시작은 /content 를 지우지 않습니다.
if subprocess.run('which conda', shell=True, capture_output=True).returncode:
    print('[중단] conda 가 없습니다. **셀 3 을 먼저 돌리십시오.**')
    print('       (노트북을 다시 업로드했다면 새 세션입니다 — 셀 3 부터입니다)')
    raise SystemExit
!conda --version
if not os.path.isdir('/content/dlbd'):
    sh('git clone -q https://github.com/nrbennet/dl_binder_design.git /content/dlbd')
sh('cd /content/dlbd && git submodule update --init --recursive 2>&1 | tail -5')
print()
print('=== include/ (환경 yml 이 여기 있습니다) ===')
sh('ls -la /content/dlbd/include/')
print()
print('=== 최상위 ===')
sh('ls /content/dlbd/')
print()
print('=== 실행 스크립트 위치 ===')
sh('find /content/dlbd -maxdepth 2 -name "*.py" | grep -Ei "interface_design|predict" ')
print()
print('=== silent_tools 위치 (버전마다 다릅니다) ===')
sh('find /content/dlbd -maxdepth 3 -type d -name silent_tools')
sh('find /content/dlbd -maxdepth 4 -name "silentfrompdbs" -o -maxdepth 4 -name "silentls" | head')
print()
print('[다음] 셀 5 가 경로와 환경 이름을 **자동으로** 찾습니다. 손으로 적지 않습니다.')""")

# ── 5 ────────────────────────────────────────────────────────────────
md("""## 5. conda 환경 두 개 만들기

셀 4 의 출력에서 본 yml 이름을 아래 두 변수에 넣으십시오.
흔한 이름은 `proteinmpnn_fastrelax.yml` 과 `af2_binder_design.yml` 입니다.

**시간이 오래 걸립니다**(각 10~30 분). PyRosetta 패키지가 1 GB 이상이고
AF2 쪽은 `tensorflow` 까지 풀어야 합니다. solve 가 멈춘 것처럼 보여도
기다리십시오.

**채널 하나를 더해야 합니다.** yml 에 적힌 채널
(pytorch / nvidia / conda-forge / defaults) 만으로는 `pyrosetta` 를 찾지
못해 44초 만에 `PackagesNotFoundError` 로 실패합니다. 아래 셀이
`conda.graylab.jhu.edu` 를 **yml 에 직접** 넣습니다 — conda 26.x 의
`conda env create` 는 `-c` 를 받지 않습니다. 실패하면 `notes/dl_binder_design_setup.md` 의 "가벼운 경로" 로
우회합니다.""")

code("""# 경로와 환경 이름을 자동으로 찾습니다.
#
# 처음에는 이것들을 네 군데에 하드코딩했고, silent_tools 가
# /content/dlbd/silent_tools 가 아니라 include/silent_tools 에 있어서
# 네 셀이 전부 실패할 상태였습니다. 추측을 흩어 놓지 않고 한 곳에서
# 탐색합니다.
import os, glob, subprocess, yaml

# --- silent_tools ---
cand = [d for d in glob.glob('/content/dlbd/**/silent_tools', recursive=True)
        if os.path.isdir(d) and glob.glob(f'{d}/silentfrompdbs')]
if not cand:
    cand = [os.path.dirname(x) for x in
            glob.glob('/content/dlbd/**/silentfrompdbs', recursive=True)]
ST = cand[0] if cand else None
print('silent_tools :', ST or '[없음] submodule 이 비었습니다')
if ST:
    print('  도구:', sorted(os.path.basename(x) for x in glob.glob(f'{ST}/silent*'))[:8])

# --- 실행 스크립트 ---
MPNN_PY = '/content/dlbd/mpnn_fr/dl_interface_design.py'
AF2_PY  = '/content/dlbd/af2_initial_guess/predict.py'
for x in (MPNN_PY, AF2_PY):
    print(('[ok] ' if os.path.isfile(x) else '[없음] ') + x)

# --- 환경 yml 과 그 안의 name: ---
print()
print('=== include/*.yml 의 환경 이름 ===')
YML = {}
for y in sorted(glob.glob('/content/dlbd/include/*.yml')):
    try:
        d = yaml.safe_load(open(y)) or {}
        nm = d.get('name', '(name 없음)')
    except Exception as e:
        nm = f'(파싱 실패: {e})'
    YML[os.path.basename(y)] = nm
    print(f'  {os.path.basename(y):32} name: {nm}')
print()
print('[확인] 위 name 이 conda 환경 이름이 됩니다. 셀 6 이 이 값을 씁니다.')""")

code("""%%time
# 두 환경을 만듭니다. 각각 수 분 ~ 15분. solve 가 멈춘 것처럼 보여도 기다리십시오.
#
# include/ 에 yml 이 셋 있을 수 있습니다 (proteinmpnn_fastrelax / af2_binder_design
# / dl_binder_design). 앞의 둘만 씁니다 — dl_binder_design.yml 은 둘을 합친
# 변형이고, 환경을 둘로 나눠야 PyRosetta(CPU)와 JAX(GPU)가 충돌하지 않습니다.
import subprocess, os, glob

# **PyRosetta 는 conda-forge 에 없습니다.** yml 의 채널 목록
# (pytorch / nvidia / conda-forge / defaults) 만으로는
#   PackagesNotFoundError: - pyrosetta
# 로 44초 만에 실패합니다. RosettaCommons 채널을 더해야 합니다.
#
# 레포가 이 채널을 yml 에 적어 두지 않은 것은 사용자가 라이선스를 직접
# 확인하도록 한 설계로 보입니다. **비상업적 사용이라는 판단은 본인의
# 것이며 독자에게 이전되지 않습니다** (레포 고지 참조).
# 채널은 **yml 에 직접** 넣습니다. conda 26.x 의 `conda env create` 는
# `-c` 를 받지 않습니다 (unrecognized arguments). yml 수정이 더 명시적이고
# 기록으로도 남습니다.
import yaml
GRAYLAB = 'https://conda.graylab.jhu.edu'
os.makedirs('/content/yml', exist_ok=True)

USE = ['proteinmpnn_fastrelax.yml', 'af2_binder_design.yml']
for y in USE:
    src = f'/content/dlbd/include/{y}'
    if not os.path.isfile(src):
        print(f'[없음] {src} — 셀 5 의 목록에서 이름을 확인하십시오'); continue
    d = yaml.safe_load(open(src))
    # conda 는 앞 채널을 우선합니다. conda-forge 보다 앞에 둡니다.
    if GRAYLAB not in d.get('channels', []):
        d['channels'] = [GRAYLAB] + d.get('channels', [])
    dst = f'/content/yml/{y}'
    # safe_dump 이 pip 하위 목록(jax[cuda12_pip]==0.4.20 등)을 보존합니다.
    yaml.safe_dump(d, open(dst, 'w'), sort_keys=False, default_flow_style=False)
    print(f'=== {y} (name: {d["name"]}) ===', flush=True)
    print('  channels:', d['channels'], flush=True)
    r = subprocess.run(f'conda env create -q -f {dst}', shell=True,
                       capture_output=True, text=True)
    print((r.stdout + r.stderr).strip()[-1500:], flush=True)

# ── AF2 환경의 JAX 조합을 **실측으로 검증된 버전**으로 맞춥니다 ──
#
# yml 의 `jax[cuda12_pip]==0.4.20` 핀은 유지되지 않습니다. 아래 줄에 버전
# 없는 `jaxlib` 이 따로 적혀 있어 pip 이 더 새 것으로 끌어올리고, 결과가
# **CPU 빌드 jaxlib** 이 됩니다:
#
#     devices [CpuDevice(id=0)]
#     An NVIDIA GPU may be present ... CUDA-enabled jaxlib is not installed
#
# 이대로 돌리면 L4 를 꽂아 둔 채 CPU 로 계산합니다 — RFdiffusion 에서
# LD_LIBRARY_PATH 때문에 27배 느려진 것과 같은 함정입니다.
#
# 2026-11-04 에 import 와 jax.devices() 로 확인한 조합:
#     jax 0.4.34 / jaxlib 0.4.34 / numpy 2.4.6 / ml_dtypes 0.6.0 / haiku 0.0.11
#
# **0.4.20 으로 내려 맞추려 하지 마십시오.** 그 시기 cuda12_pip wheel 이
# find-links 에서 사라졌고, numpy 를 1.x 로 내리면 이번엔 tensorflow 2.21 이
# ml_dtypes>=0.5.1 을 요구해 또 충돌합니다. 실제로 그 경로를 한 번 갔고
# np.trapz -> take_current_trace -> ml_dtypes 로 세 번 연쇄했습니다.
AF2_ENV_NAME  = YML.get('af2_binder_design.yml', 'af2_binder_design')
MPNN_ENV_NAME = YML.get('proteinmpnn_fastrelax.yml', 'proteinmpnn_binder_design')

# **numpy 는 두 환경에서 서로 반대로 고정해야 합니다.**
#   AF2  : tensorflow 2.21 / ml_dtypes>=0.5.1 이 numpy>=2 를 요구
#   MPNN : torch 2.0 이 numpy 1.x ABI 로 빌드됨
#          -> numpy 2 면 "Failed to initialize NumPy: _ARRAY_API not found"
# 환경을 둘로 나눈 설계가 여기서 값을 합니다. 한 환경이면 풀 수 없습니다.
print()
print('=== MPNN 환경: numpy 를 1.x 로 (torch 2.0 ABI) ===')
r = subprocess.run(f'conda run -n {MPNN_ENV_NAME} pip install -q "numpy<2"',
                   shell=True, capture_output=True, text=True)
print((r.stdout + r.stderr).strip()[-600:])

print()
print('=== AF2 환경의 JAX 를 CUDA 빌드로 맞춥니다 (wheel 2GB+, 5~10분) ===')
r = subprocess.run(
    f'conda run -n {AF2_ENV_NAME} pip install -q '
    f'"jax[cuda12]==0.4.34" "jaxlib==0.4.34" '
    f'"ml-dtypes>=0.5.1" "numpy>=2" "dm-haiku==0.0.11"',
    shell=True, capture_output=True, text=True)
print((r.stdout + r.stderr).strip()[-1200:])

# **AlphaFold 소스 한 줄 패치.** mmcif_parsing.py 가 Bio.Data.SCOPData 를
# import 하는데 biopython 1.83 에는 없습니다(PDBData 로 이동). 버전을 내려
# 맞추려 하지 마십시오 — 1.79 는 Python 3.11 빌드가 없습니다.
# **이것은 레포 소스를 수정하는 것입니다.** 재현하는 사람이 알아야 합니다.
MM = '/content/dlbd/af2_initial_guess/alphafold/data/mmcif_parsing.py'
src = open(MM).read()
if 'PDBData' not in src:
    # 줄바꿈은 chr(10) 으로 만듭니다 — 생성기의 삼중인용 문자열 안에서
    # 역슬래시-n 을 쓰면 **생성 시점에** 실제 줄바꿈이 되어 노트북 셀의
    # 문자열이 깨집니다 (실제로 깨져 있었습니다).
    SHIM = [
        'try:',
        '  from Bio.Data import SCOPData  # biopython < 1.83',
        'except ImportError:',
        '  from Bio.Data import PDBData as _P',
        '  class SCOPData:  # shim',
        '    protein_letters_3to1 = _P.protein_letters_3to1_extended',
    ]
    open(MM, 'w').write(src.replace('from Bio.Data import SCOPData',
                                    chr(10).join(SHIM), 1))
    print()
    print('[패치] mmcif_parsing.py — SCOPData -> PDBData shim')

print()
# ProteinMPNN 은 이 버전에서 서브모듈이 아닙니다(.gitmodules 가 없습니다).
# importtest 가 이것을 import 하므로 **시험 전에** 받아야 합니다.
MPNN_DIR = os.path.join(os.path.dirname(MPNN_PY), 'ProteinMPNN')
if not os.path.isdir(MPNN_DIR):
    print('ProteinMPNN 이 없습니다. clone 합니다.')
    subprocess.run(f'cd {os.path.dirname(MPNN_PY)} && '
                   'git clone -q https://github.com/dauparas/ProteinMPNN.git',
                   shell=True)
print('ProteinMPNN .pt:',
      len(glob.glob(f'{MPNN_DIR}/**/*.pt', recursive=True)), '개')
print()
print('=== 레포가 제공하는 import 시험 (include/importtests) ===')
print('    README: "여기서 에러가 나면 설치가 잘못된 것"')
for t, e in (('proteinmpnn_importtest.py', MPNN_ENV_NAME),
             ('af2_importtest.py', AF2_ENV_NAME)):
    r = subprocess.run(
        f'cd /content/dlbd && conda run -n {e} python include/importtests/{t}',
        shell=True, capture_output=True, text=True)
    hits = [l for l in (r.stdout + r.stderr).splitlines()
            if any(k in l.lower() for k in ('pass', 'fail', 'gpu', 'error'))]
    print(f"[{'ok' if r.returncode == 0 else '실패'}] {t} ({e})")
    print('   ' + (chr(10) + '   ').join(hits[-4:]))

print()
print('=== 직접 import 확인 ===')
for env, code in (
    (MPNN_ENV_NAME,
     "import numpy, pyrosetta, torch, Bio;"
     " print('  numpy', numpy.__version__, '/ torch', torch.__version__);"
     " print('  torch<->numpy:', torch.from_numpy(numpy.zeros(3)).sum().item())"),
    (AF2_ENV_NAME,
     "import numpy, ml_dtypes, jax, haiku, pyrosetta, Bio;"
     " print('  numpy', numpy.__version__, '/ ml_dtypes', ml_dtypes.__version__);"
     " print('  jax', jax.__version__, '/ haiku', haiku.__version__);"
     " print('  jax devices:', jax.devices())"),
):
    r = subprocess.run(f'conda run -n {env} python -c "{code}"',
                       shell=True, capture_output=True, text=True)
    print(f"[{'ok' if r.returncode == 0 else '실패'}] {env}")
    print((r.stdout + r.stderr).strip()[-900:])
print()
print('[확인] jax devices 가 [CudaDevice(id=0)] 이어야 합니다.')
print('       CpuDevice 만 나오면 **여기서 멈추십시오** — AF2 가 조용히')
print('       CPU 로 돌아 100건이 몇십 시간이 됩니다.')
print()
print('[남은 위험] 레포의 vendored AlphaFold 코드가 jax 0.4.20 시절 API 를')
print('       쓰면 0.4.34 에서 깨질 수 있습니다. 설치로는 알 수 없고')
print('       셀 8(게이트 B)이 실제 시험입니다.')

print()
!conda env list
print()
print('[확인] 위 목록에 환경 두 개가 보여야 합니다. 그 이름을 셀 7 이 씁니다.')""")

md("""### 환경 이름 고정

`conda env list` 에서 본 이름을 아래가 자동으로 집습니다. 뒤 셀들은 모두
이 두 변수를 씁니다 — 이름을 여러 곳에 적지 않습니다.""")

code("""import subprocess, re

r = subprocess.run('conda env list', shell=True, capture_output=True, text=True)
names = [l.split()[0] for l in r.stdout.splitlines()
         if l.strip() and not l.startswith('#')]
print('환경:', names)

def pick(pats):
    for p in pats:
        for n in names:
            if p in n.lower():
                return n
    return None

MPNN_ENV = pick(['mpnn', 'fastrelax', 'proteinmpnn'])
AF2_ENV  = pick(['af2', 'alphafold'])
print()
print('MPNN_ENV =', MPNN_ENV)
print('AF2_ENV  =', AF2_ENV)
if not (MPNN_ENV and AF2_ENV):
    print()
    print('[중단] 환경을 찾지 못했습니다. 위 목록을 보고 아래처럼 직접 넣으십시오:')
    print("  MPNN_ENV = 'proteinmpnn_binder_design'")
    print("  AF2_ENV  = 'af2_binder_design'")
else:
    print()
    print('[ok] 다음 셀로')""")

# ── 6 ────────────────────────────────────────────────────────────────
md("""## 6. 가중치 받기 (ProteinMPNN + AF2)

**ProteinMPNN 가중치는 서브모듈로 오지 않습니다**(이 버전에는 `.gitmodules`
자체가 없습니다). `mpnn_fr/` 안에 직접 clone 해야 합니다 — 레포 README 가
지시하는 방식입니다. AF2 파라미터는 따로 받습니다(약 3.5 GB);
`af2_initial_guess` 는 `model_1_ptm` 하나만 씁니다.

이 셀은 **`BASE` 환경변수 묶음도 정의합니다.** 뒤의 모든 셀이 이것을 씁니다.
두 가지가 들어 있고 둘 다 실패를 겪고 추가한 것입니다 —
`silent_tools` 의 셰뱅이 `#!/usr/bin/env python` 이라 conda base 를 집어
`No module named numpy` 가 나므로 **환경의 bin 을 PATH 앞에** 붙이고,
`PYTHONPATH` 는 `silent_tools` 의 **부모 디렉터리**를 가리켜야 합니다
(디렉터리 자체를 주면 동명의 `silent_tools.py` 가 패키지를 가립니다).""")

code("""%%time
import os, glob, subprocess

def sh(c):
    r = subprocess.run(c, shell=True, capture_output=True, text=True)
    print((r.stdout + r.stderr).strip()[-1500:])

# --- 뒤의 모든 셀이 쓰는 환경 ---
PATH0 = f'/usr/local/envs/{MPNN_ENV}/bin:{ST}:' + os.environ['PATH']
BASE = dict(os.environ, PATH=PATH0, PYTHONPATH=os.path.dirname(ST))
print('BASE PATH 선두  :', BASE['PATH'].split(':')[0])
print('BASE PYTHONPATH :', BASE['PYTHONPATH'])
w = subprocess.run('which python', shell=True, capture_output=True,
                   text=True, env=BASE).stdout.strip()
print('which python    :', w)
print()

print('=== ProteinMPNN 가중치 ===')
MPNN_DIR = os.path.join(os.path.dirname(MPNN_PY), 'ProteinMPNN')
if not os.path.isdir(MPNN_DIR):
    print('없습니다. clone 합니다 (레포 README 가 지시하는 방식).')
    sh(f'cd {os.path.dirname(MPNN_PY)} && '
       'git clone -q https://github.com/dauparas/ProteinMPNN.git')
pt = glob.glob(f'{MPNN_DIR}/**/*.pt', recursive=True)
print(f'.pt {len(pt)}개')
have = [x for x in pt if x.endswith('v_48_020.pt')]
print('[ok] v_48_020.pt' if have else '[실패] v_48_020.pt 가 없습니다')

PARAMS = '/content/dlbd/af2_initial_guess/model_weights/params'
need = f'{PARAMS}/params_model_1_ptm.npz'
print()
print('=== AF2 파라미터 ===')
if os.path.isfile(need):
    print('[ok] 이미 있습니다:', need)
else:
    os.makedirs(PARAMS, exist_ok=True)
    sh(f'cd {PARAMS} && wget -q --show-progress '
       'https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar '
       '&& tar --extract --file=alphafold_params_2022-12-06.tar '
       '&& rm alphafold_params_2022-12-06.tar')
    print('[ok]' if os.path.isfile(need) else '[실패]', need)""")
md("""## 7. silent 파일 만들기

PDB 묶음을 silent 하나로 묶습니다. **설계군·대조군·게이트를 따로** 만듭니다
— 섞으면 어느 설계가 어느 군인지 뒤에서 가릴 수 없습니다.""")

code("""import os, glob, shutil, subprocess
os.makedirs('/content/h1', exist_ok=True)
os.chdir('/content/h1')

ROOT = '/content/drive/MyDrive/RFdiffusion'
GROUPS = {
    'epi':  sorted(glob.glob(f'{ROOT}/outputs/pdl1_epi_*.pdb')) +
            sorted(glob.glob(f'{ROOT}/outputs/pdl1_epi2_*.pdb')),
    'ctrl': sorted(glob.glob(f'{ROOT}/outputs/pdl1_ctrl_*.pdb')) +
            sorted(glob.glob(f'{ROOT}/outputs/pdl1_ctrl2_*.pdb')),
    'gate': [f'{ROOT}/inputs/PD1_PDL1_binderfirst.pdb'],
}

env = dict(os.environ, PATH=f"{ST}:{os.environ['PATH']}")   # ST 는 셀 5 에서 찾았습니다

for name, files in GROUPS.items():
    d = f'/content/h1/pdb_{name}'
    os.makedirs(d, exist_ok=True)
    for f in files:
        shutil.copy(f, d)
    r = subprocess.run(
        f'cd {d} && {ST}/silentfrompdbs *.pdb > /content/h1/{name}.silent',
        shell=True, capture_output=True, text=True, env=env)
    n = os.path.getsize(f'/content/h1/{name}.silent') if \\
        os.path.isfile(f'/content/h1/{name}.silent') else 0
    print(f'{name:5} PDB {len(files):3}개 -> silent {n/1e6:6.2f} MB'
          + ('' if n else f'   [실패] {(r.stdout + r.stderr)[-300:]}'))

print()
!ls -la /content/h1/*.silent
print()
print('[확인] 세 silent 가 모두 0 보다 커야 합니다.')""")


code("""# silent 의 태그를 셉니다. 바이트 크기만 보면 조용한 누락을 놓칩니다.
import subprocess, os

def tags_of(silent):
    # stdout 만 읽습니다 (stderr 를 섞으면 경고 단어가 태그로 세어집니다).
    r = subprocess.run(f'{ST}/silentls {silent}', shell=True,
                       capture_output=True, text=True, env=BASE)
    return [t for t in r.stdout.split() if t]

for name, want in (('epi', 50), ('ctrl', 50), ('gate', 1)):
    t = tags_of(f'/content/h1/{name}.silent')
    print(f'{name:5} 태그 {len(t):3} / 기대 {want:3}'
          + ('  [ok]' if len(t) == want else '  [멈추십시오]'))""")
# ── 8 ────────────────────────────────────────────────────────────────
md("""## 8. 게이트 B — native PD-1 로 `pae_interaction` 이 작동하는가

**설계군을 돌리기 전에 이것부터.** 실제로 결합하는 복합체에서
`pae_interaction` 이 10 보다 크면 이 지표로 설계를 판정할 수 없습니다.

여기서는 ProteinMPNN 을 거치지 않습니다 — native 서열을 그대로 AF2 에
넣어 재예측합니다. 사슬 A = PD-1(binder), B = PD-L1(타깃) 입니다.""")

code("""%%time
import subprocess, os, glob

print('predict.py:', os.path.isfile(AF2_PY))

cmd = (f'conda run -n {AF2_ENV} python {AF2_PY} '
       f'-silent /content/h1/gate.silent '
       f'-outsilent /content/h1/gate_af2.silent '
       f'-scorefilename /content/h1/gate_af2.sc 2>&1 | tail -30')
print('$', cmd)
r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
print((r.stdout + r.stderr).strip()[-3000:])""")

code("""# 게이트 판정
import os, re

sc = '/content/h1/gate_af2.sc'
if not os.path.isfile(sc):
    print('[중단] 점수 파일이 없습니다. 위 셀의 오류를 먼저 해결하십시오.')
else:
    lines = [l.split() for l in open(sc) if l.strip()]
    hdr = lines[0]
    print('열:', hdr)
    print()
    try:
        i = hdr.index('pae_interaction')
    except ValueError:
        i = next((k for k, h in enumerate(hdr) if 'pae' in h.lower()), None)
        print(f'[주의] pae_interaction 열을 찾지 못해 {hdr[i]!r} 를 씁니다')
    vals = []
    for row in lines[1:]:
        if len(row) > i:
            try:
                vals.append((float(row[i]), row[-1]))
            except ValueError:
                pass
    for v, n in vals:
        print(f'  {n}: pae_interaction = {v:.2f}')
    print()
    if vals and min(v for v, _ in vals) < 10:
        print('  [통과] native 복합체가 10 미만입니다. 설계군을 돌릴 수 있습니다.')
        print('  주의: PD-1 사슬의 84->93 결손 8잔기(85-92)를 메우지 않은 값입니다.')
        print('  사전 등록 §4 에 따라 메운 값도 함께 보고해야 합니다.')
    else:
        print('  [실패] native 복합체가 10 이상입니다. **여기서 멈추십시오.**')
        print('  native 가 통과하지 못하는 기준으로 설계를 판정할 수 없습니다.')
        print('  H1 을 dG_separated 만으로 축소하고 그 사실을 README 한계에')
        print('  적으십시오. (Rosetta 게이트 A 는 2026-11-04 통과: -32.74 REU)')""")

# ── 9 ────────────────────────────────────────────────────────────────
md("""## 9. ProteinMPNN — binder 서열 설계

백본당 서열 **1 개**입니다(`-seqs_per_struct 1`). 군 간 비교에는 설계당
하나가 충분하고, 모자라면 늘리면 됩니다. 8 개로 하면 AF2 시간이 8 배입니다.

타깃 사슬은 고정됩니다 — 설계하는 것은 binder 사슬뿐입니다.
이 단계는 **CPU** 이고 FastRelax 때문에 느립니다. 먼저 시간을 재십시오.""")

code("""%%time
import subprocess, os

print('dl_interface_design.py:', os.path.isfile(MPNN_PY))

# 먼저 2건으로 비용을 잽니다. 100건을 눈 감고 돌리지 않습니다.
!{'conda run -n ' + MPNN_ENV + ' python ' + MPNN_PY + ' -h 2>&1 | head -30'}""")

code("""%%time
# 비용 측정: 설계군에서 2건만
import subprocess, time, os


os.makedirs('/content/h1/probe', exist_ok=True)
r = subprocess.run(
    f'{ST}/silentslice /content/h1/epi.silent $({ST}/silentls '
    f'/content/h1/epi.silent | head -2) > /content/h1/probe/epi2.silent',
    shell=True, capture_output=True, text=True,
    env=dict(os.environ, PATH=f"{ST}:{os.environ['PATH']}"))
print((r.stdout + r.stderr)[-500:])
print('probe silent:', os.path.getsize('/content/h1/probe/epi2.silent'))

t = time.time()
cmd = (f'conda run -n {MPNN_ENV} python {MPNN_PY} '
       f'-silent /content/h1/probe/epi2.silent '
       f'-outsilent /content/h1/probe/epi2_mpnn.silent '
       f'-seqs_per_struct 1 2>&1 | tail -20')
r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
print((r.stdout + r.stderr)[-2000:])
dt = time.time() - t
print()
print(f'2건에 {dt/60:.1f}분 -> 100건이면 약 {dt/2*100/60:.1f}분')""")

code("""%%time
# 본실행: **병렬** MPNN. FastRelax 가 단일 스레드이므로 CPU 를 쪼개면
# 5 배 빨라집니다 (실측 134분 -> 25.7분). GPU 는 오히려 비웁니다 —
# 6 프로세스가 한 GPU 를 다투면 느려집니다.
import subprocess, os, glob, math
from concurrent.futures import ThreadPoolExecutor

WORK = '/content/h1/mpnn'
os.makedirs(WORK, exist_ok=True)
NCHUNK = NWORKER = 6

def tags_of(silent):
    # stdout 만 읽습니다. stderr 를 섞으면 silentls 의 경고 단어가
    # 태그로 세어집니다 (실제로 56/57 로 어긋났습니다).
    r = subprocess.run(f'{ST}/silentls {silent}', shell=True,
                       capture_output=True, text=True, env=BASE)
    return [t for t in r.stdout.split() if t]

def merge_silent(parts, out):
    # 2 번째 조각부터 파일 헤더 두 줄을 뗍니다.
    # silentcat 같은 명령은 **존재하지 않습니다**.
    with open(out, 'w') as o:
        for k, p in enumerate(parts):
            for ln in open(p):
                if k and (ln.startswith('SEQUENCE:') or
                          (ln.startswith('SCORE:') and
                           ln.rstrip().split()[-1] == 'description')):
                    continue
                o.write(ln)

def run_chunk(job):
    g, i, tags = job
    sl = f'{WORK}/{g}_{i}.silent'
    out = f'{WORK}/{g}_{i}_mpnn.silent'
    joined = ' '.join(tags)
    subprocess.run(f'{ST}/silentslice /content/h1/{g}.silent {joined} > {sl}',
                   shell=True, env=BASE)
    r = subprocess.run(
        f'conda run -n {MPNN_ENV} python {MPNN_PY} -silent {sl} '
        f'-outsilent {out} -seqs_per_struct 1 '
        f'-checkpoint_name {WORK}/{g}_{i}.point',
        shell=True, capture_output=True, text=True,
        env=dict(BASE, CUDA_VISIBLE_DEVICES=''))
    return g, i, os.path.isfile(out), (r.stdout + r.stderr)[-300:]

jobs = []
for g in ('epi', 'ctrl'):
    t = tags_of(f'/content/h1/{g}.silent')
    k = math.ceil(len(t) / NCHUNK)
    nch = math.ceil(len(t) / k)
    jobs += [(g, i, t[i * k:(i + 1) * k]) for i in range(nch)]
    print(f'{g}: 입력 태그 {len(t)} -> {nch}조각')

with ThreadPoolExecutor(NWORKER) as ex:
    for g, i, ok, log in ex.map(run_chunk, jobs):
        print(f'  {g}_{i}: ' + ('ok' if ok else '실패'))
        if not ok:
            print('   ', log)""")

code("""%%time
# 누락 검증 -> 재실행 -> 합치기. **건너뛰지 마십시오.**
# MPNN 은 오류를 내지 않고 설계를 떨어뜨립니다 (실측 100건 중 3건).
# 세지 않으면 48 대 49 로 AF2 를 2시간 돌리고도 n 이 깨진 것을 모릅니다.
import subprocess, os, glob

for g in ('epi', 'ctrl'):
    want = tags_of(f'/content/h1/{g}.silent')
    parts = sorted(glob.glob(f'{WORK}/{g}_*_mpnn.silent'))
    have = {t.split('_dldesign')[0] for p in parts for t in tags_of(p)}
    missing = [t for t in want if t not in have]
    print(f'{g}: 입력 {len(want)}  출력 {len(have)}  '
          f'누락 {len(missing)} {missing}')
    if missing:
        d = f'{WORK}/retry_{g}'
        os.makedirs(d, exist_ok=True)
        joined = ' '.join(missing)
        subprocess.run(f'{ST}/silentslice /content/h1/{g}.silent '
                       f'{joined} > {d}/in.silent', shell=True, env=BASE)
        r = subprocess.run(
            f'conda run -n {MPNN_ENV} python {MPNN_PY} -silent {d}/in.silent '
            f'-outsilent {d}/out.silent -seqs_per_struct 1 -debug '
            f'-checkpoint_name {d}/retry.point',
            shell=True, capture_output=True, text=True,
            env=dict(BASE, CUDA_VISIBLE_DEVICES=''))
        print((r.stdout + r.stderr)[-600:])
        if os.path.isfile(f'{d}/out.silent'):
            parts.append(f'{d}/out.silent')
    merged = f'/content/h1/{g}_mpnn.silent'
    merge_silent(parts, merged)
    n = len(set(tags_of(merged)))
    verdict = '[ok]' if n == len(want) else '[n 이 깨졌습니다 — 멈추십시오]'
    print(f'  -> {merged}  고유 태그 {n}  {verdict}')""")
# ── 10 ───────────────────────────────────────────────────────────────
md("""## 10. AF2 initial guess — 복합체 재예측

여기가 H1 의 핵심입니다. ProteinMPNN 이 올린 서열로 복합체를 **다른
모델이** 다시 예측하고, 그 예측이 RFdiffusion 의 백본과 맞는지 봅니다.
맞지 않으면 RFdiffusion 의 주장이 자기만의 환상이었다는 뜻입니다.

GPU 단계입니다. 먼저 2건으로 시간을 재십시오.""")

code("""%%time
# 비용 측정: 2건
import subprocess, time, os
env = dict(os.environ, PATH=f"{ST}:{os.environ['PATH']}")

subprocess.run(
    f'{ST}/silentslice /content/h1/epi_mpnn.silent $({ST}/silentls '
    f'/content/h1/epi_mpnn.silent | head -2) > /content/h1/probe/epi2_m.silent',
    shell=True, env=env)

t = time.time()
r = subprocess.run(
    f'conda run -n {AF2_ENV} python {AF2_PY} '
    f'-silent /content/h1/probe/epi2_m.silent '
    f'-outsilent /content/h1/probe/epi2_af2.silent '
    f'-scorefilename /content/h1/probe/epi2_af2.sc 2>&1 | tail -20',
    shell=True, capture_output=True, text=True)
print((r.stdout + r.stderr)[-2000:])
dt = time.time() - t
print()
print(f'2건에 {dt/60:.1f}분 -> 100건이면 약 {dt/2*100/60:.1f}분')
print('(첫 실행에는 JAX 컴파일 시간이 포함되므로 100건 추정은 과대평가입니다)')""")

code("""%%time
# 본실행: 군별로 따로 돕니다. 출력·점수·체크포인트를 **Drive 에** 두어
# 세션이 끊겨도 잃지 않습니다 — 같은 셀을 다시 돌리면 체크포인트가 끝난
# 태그를 건너뛰고 이어갑니다 (여기서는 그 건너뛰기가 원하는 동작입니다).
import subprocess, os, time

DEST = '/content/drive/MyDrive/RFdiffusion/h1'
os.makedirs(DEST, exist_ok=True)

for g in ('epi', 'ctrl'):
    print(f'=== {g} (실측 50건 59분, L4) ===', flush=True)
    t0 = time.time()
    r = subprocess.run(
        f'cd {os.path.dirname(AF2_PY)} && conda run -n {AF2_ENV} python {AF2_PY} '
        f'-silent /content/h1/{g}_mpnn.silent '
        f'-outsilent {DEST}/{g}_af2.silent '
        f'-scorefilename {DEST}/{g}_af2.sc '
        f'-checkpoint_name {DEST}/{g}_af2.point',
        shell=True, capture_output=True, text=True, env=BASE)
    print((r.stdout + r.stderr)[-800:])
    print(f'소요 {(time.time() - t0) / 60:.1f}분', flush=True)
    sc = f'{DEST}/{g}_af2.sc'
    if os.path.isfile(sc):
        n = sum(1 for l in open(sc) if l.startswith('SCORE:')) - 1
        print(f'점수 {n}건' + ('' if n == 50 else '   [확인 필요]'))""")
# ── 11 ───────────────────────────────────────────────────────────────
md("""## 11. H1 판정 + PDB 추출 + Drive 백업

**2026-11-04 변경 — `pae_interaction` 으로 선별하지 않습니다.**
이 지표의 양성 대조군(native PD-1)이 26.06 으로 실패했으므로 사전 등록 §4 에
따라 "측정 실패"로 보고하고 순위화에 쓰지 않습니다. 선별에 쓰면 쓰지 않기로
판정한 지표가 **선택 편향으로 되돌아옵니다.** 100 건을 전부 Rosetta 로
넘기고, 판정은 승격된 `dG_separated` 가 내립니다.""")

code("""import os, statistics as st

def read_sc(path):
    # 헤더 이름으로 파싱합니다. 열 순서를 가정하지 않습니다.
    rows, hdr = [], None
    for ln in open(path):
        if not ln.startswith('SCORE:'):
            continue
        f = ln.split()[1:]
        if hdr is None:
            hdr = f
            continue
        rows.append(dict(zip(hdr, f)))
    return rows

DEST = '/content/drive/MyDrive/RFdiffusion/h1'
for g, label in (('epi', '설계군'), ('ctrl', '대조군')):
    p = f'{DEST}/{g}_af2.sc'
    if not os.path.isfile(p):
        print(f'[{g}] 점수 파일 없음')
        continue
    rows = read_sc(p)
    pae = sorted(float(r['pae_interaction']) for r in rows)
    print(f'== {label}  n={len(rows)} ==')
    print(f'  pae_interaction  중앙값 {st.median(pae):.2f}'
          f'   최저 {pae[0]:.2f}   최고 {pae[-1]:.2f}')
    print(f'  <10 {sum(x < 10 for x in pae)}건   '
          f'<15 {sum(x < 15 for x in pae)}건   '
          f'<20 {sum(x < 20 for x in pae)}건')
    for r in sorted(rows, key=lambda r: float(r['pae_interaction']))[:5]:
        tag = r['description'][:28]
        print(f'    {tag:30s} pae_int {float(r["pae_interaction"]):6.2f}'
              f'  pLDDT_b {float(r["plddt_binder"]):5.1f}'
              f'  rmsd {float(r["binder_aligned_rmsd"]):5.2f}')
    print()

print('이 숫자는 판정 근거가 아닙니다 (사전 등록 §4). 관측값입니다.')""")

code("""# 100 건 전부 PDB 로 꺼내 Drive 에 백업
import os, subprocess, glob

OUT = '/content/h1/af2_pdb'
DEST = '/content/drive/MyDrive/RFdiffusion/h1'
for g in ('epi', 'ctrl'):
    d = f'{OUT}/{g}'
    os.makedirs(d, exist_ok=True)
    r = subprocess.run(f'cd {d} && {ST}/silentextract {DEST}/{g}_af2.silent',
                       shell=True, capture_output=True, text=True, env=BASE)
    pdbs = glob.glob(f'{d}/*.pdb')
    ok = '' if len(pdbs) == 50 else '   [확인 필요]'
    print(f'{g}: PDB {len(pdbs)}개{ok}')
    if len(pdbs) != 50:
        print((r.stdout + r.stderr)[-500:])
    tar = f'{DEST}/{g}_af2_pdb.tar.gz'
    subprocess.run(f'cd {OUT} && tar czf {tar} {g}', shell=True)
    print(f'   -> {tar}  {os.path.getsize(tar) / 1e6:.1f} MB')

print()
print('[다음] 두 tar.gz 를 내려받아 outputs/h1_af2/ 에 풀고, WSL2 에서')
print('scripts/h1_rosetta.py 를 군별로 돌린 뒤 --compare 합니다.')
print('python 에 -u 를 주십시오 — 빼면 로그가 8 KB 단위로 버퍼링되어')
print('한 시간 넘게 진행을 볼 수 없습니다 (실제로 겪음).')
print()
print('  게이트 A 는 2026-11-04 통과했습니다 (dG_separated -32.74 REU).')
print('  실측 1.8분/건 -> 100건에 약 3시간.')""")

# ──────────────────────────────────────────────────────────────────────
nb = {
    "cells": [
        {"cell_type": t,
         "metadata": {},
         "source": (src.rstrip("\n") + "\n").splitlines(keepends=True),
         **({"outputs": [], "execution_count": None} if t == CODE else {})}
        for t, src in CELLS
    ],
    "metadata": {
        "accelerator": "GPU",
        "colab": {"provenance": [], "toc_visible": True},
        "kernelspec": {"display_name": "Python 3", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 0,
}

dst = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "H1_dl_binder_design.ipynb")
with io.open(dst, "w", encoding="utf-8", newline="\n") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
    f.write("\n")
print(f"{dst}  ({len(CELLS)}셀: "
      f"markdown {sum(t == MD for t, _ in CELLS)}, "
      f"code {sum(t == CODE for t, _ in CELLS)})")
