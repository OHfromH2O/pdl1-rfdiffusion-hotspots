# Colab에서 RFdiffusion이 조용히 CPU로 도는 문제

**2026-11-03 확인.** 커뮤니티 RFdiffusion Colab 노트북에서 발생합니다.
경고도 에러도 없이 **27배 느려집니다.** 측정하지 않으면 알 수 없습니다.

---

## 증상

| | |
|---|---|
| 겉보기 | 정상 실행. 결과 PDB도 제대로 나옴 |
| 속도 | 59초/스텝 (205잔기, T=25 → 24.7분/설계) |
| 정상 속도 | **2.2초/스텝** (→ 약 2분/설계, T=50) |
| `nvidia-smi` | 실행 내내 `0 %, 3 MiB` |

T=50 × 100 설계 기준 **90시간 vs 3.5시간**입니다. 무료 T4로 가능한 작업이
불가능한 작업으로 보이게 됩니다.

---

## 원인

노트북이 `run_inference.py`를 서브프로세스로 띄울 때 환경을 이렇게 만듭니다.

```python
SP = f"{ENV}/lib/python3.10/site-packages"
LDPATH = ":".join(sorted(glob.glob(f"{SP}/nvidia/*/lib")) + [f"{SP}/torch/lib"])
...
proc = subprocess.Popen(
    cmd, cwd="/content/RFdiffusion",
    env=dict(os.environ, DGLBACKEND="pytorch", LD_LIBRARY_PATH=LDPATH), ...)
```

`LDPATH`에 **`/usr/lib64-nvidia`가 없습니다.** 그리고
`dict(os.environ, LD_LIBRARY_PATH=LDPATH)`는 상속값을 **덮어씁니다.**

Colab에서 `libcuda.so.1`(드라이버 라이브러리)은 `/usr/lib64-nvidia`에만 있고,
`ldconfig` 캐시에는 등재되어 있지 않습니다. 커널의 `LD_LIBRARY_PATH`가
그 경로 **하나뿐**이라 덮어쓰면 드라이버가 사라집니다.

그러면 `torch.cuda.is_available()`이 `False`를 돌려주고,
RFdiffusion은 설계대로 CPU를 고릅니다.

```python
# inference/model_runners.py
if torch.cuda.is_available():
    self.device = torch.device('cuda')
else:
    self.device = torch.device('cpu')   # ← 여기로 떨어진다
```

**양쪽 경로가 모두 필요합니다.**

| 라이브러리 | 어디에 | 누가 쓰는가 |
|---|---|---|
| `libcuda.so.1` | `/usr/lib64-nvidia` | torch (드라이버 접근) |
| `libcusparse.so.11` 등 | `{SP}/nvidia/*/lib` | dgl 1.1.3+cu118 |

하나만 넣으면 각각 다른 방식으로 실패합니다 — 드라이버 경로만 넣으면
`run_inference.py`가 `OSError: libcusparse.so.11`로 **죽고**, pip 경로만 넣으면
**죽지 않고 CPU로 돕니다.** 후자가 위험합니다.

---

## 수정

`subprocess.Popen` 호출 바로 앞에 한 줄을 넣습니다.

```python
LDPATH = LDPATH + ":" + ":".join(
    p for p in os.environ.get("LD_LIBRARY_PATH", "").split(":") if p)
```

상속값을 뒤에 덧붙입니다. pip 경로가 앞이라 우선순위는 그대로 유지됩니다.

**주의: 이 수정은 노트북 셀 안에만 존재합니다.** 런타임이 끊기면 사라지고,
다음 실행은 다시 조용히 CPU로 돕니다. 세션마다 확인해야 합니다.

---

## 확인 방법

`run_inference.py` 최상단에 넣어 두면 매 실행의 첫 줄에 상태가 찍힙니다.

```python
import os, torch
print('### TOP cuda:', torch.cuda.is_available(),
      '| CVD=', repr(os.environ.get('CUDA_VISIBLE_DEVICES')),
      '| torch.version.cuda=', torch.version.cuda, flush=True)
try:
    torch.cuda.init(); print('### TOP init ok, count=', torch.cuda.device_count(), flush=True)
except Exception as e:
    print('### TOP init ERR:', type(e).__name__, e, flush=True)
```

고장 상태:
```
### TOP cuda: False | CVD= None | torch.version.cuda= 11.8
### TOP init ERR: RuntimeError Found no NVIDIA driver on your system.
```

정상 상태:
```
### TOP cuda: True | CVD= None | torch.version.cuda= 11.8
### TOP init ok, count= 1
```

---

## 진단 과정에서 배운 것

이 문제를 찾는 데 약 2.5시간이 들었습니다. 헛디딘 지점을 남겨 둡니다.

**1. `python -c` 단독 테스트는 이 버그를 못 잡습니다.**
`/content/rfenv/bin/python`을 그냥 실행하면 커널의 `LD_LIBRARY_PATH`를
상속해 `cuda: True`가 나옵니다. **노트북이 실제로 만드는 환경을 그대로
재현해야** 합니다 — 같은 `env=` 딕셔너리로 띄워야 보입니다.

**2. `subprocess`의 stdout은 Colab 셀 출력에 안 올라옵니다.**
`capture_output=True, text=True` 후 직접 `print`해야 합니다.
그러지 않으면 `CompletedProcess(...)` repr만 보고 아무것도 알 수 없습니다.

**3. `nvidia-smi` 로그를 `sort | uniq -c`로 뭉개면 안 됩니다.**
실행 전후의 유휴 구간이 섞여 들어가 `0 %`가 과대 집계됩니다.
처음 이 집계를 보고 CPU 폴백이라 판단했는데, 그 샘플은 대부분
건너뛴 실행 뒤의 유휴 시간이었습니다. **시계열 그대로 보십시오.**

**4. 모니터 자체를 먼저 검증해야 합니다.**
확실히 GPU를 쓰는 작업(큰 행렬곱 25초)을 돌려 `nvidia-smi`가 `100 %`를
보는지 확인한 뒤에야 `0 %`가 의미를 갖습니다.

**5. 출력 이름을 재사용하면 cautious mode가 조용히 건너뜁니다.**
```
(cautious mode) Skipping this design because /content/outputs/X_0.pdb already exists.
```
0.8분에 끝난 "실행"의 GPU 측정값을 한참 들여다봤습니다.
측정할 때는 **매번 새 이름**을 쓰십시오.

**6. 노트북의 Fast debug 모드는 쓸 수 없습니다.**
`diffuser.T=5`를 넘기는데 RFdiffusion이 `T >= 15`를 단정합니다.
```
AssertionError: With discrete time and T < 15, the schedule is badly approximated
```
빠른 확인은 Export check(크기를 줄임) 쪽을 쓰거나 `iterations=25`로 하십시오.

---

## 측정 환경

| | |
|---|---|
| 플랫폼 | Google Colab 무료 티어 |
| GPU | Tesla T4, 드라이버 580.82.07 |
| rfenv | Python 3.10, torch 2.0.1+cu118, dgl 1.1.3+cu118 |
| RFdiffusion | Complex_base_ckpt.pt |
| 입력 | PD-L1 4ZQK 사슬 A (115잔기) + 설계 70-100잔기 ≈ 205잔기 |
| 측정일 | 2026-11-03 |

**참고 — T4의 정상 속도.** 2.2초/스텝은 205잔기에서 타당한 값입니다.
A100이 100잔기에서 약 3초/스텝으로 보고되고, attention이 잔기 수의 2차이므로
205잔기는 약 3~4배, T4는 A100의 1/8~1/12입니다. 만약 측정값이 60초 근처라면
**그 자체로 CPU를 의심할 근거**입니다.
