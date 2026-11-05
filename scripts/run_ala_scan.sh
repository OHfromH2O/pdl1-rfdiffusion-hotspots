#!/usr/bin/env bash
# WSL2에서 Ala 스캔을 돌리기 위한 래퍼.
#
# 왜 래퍼가 필요한가
# ------------------
# PyRosetta 는 네이티브 Windows 에서 동작하지 않습니다. 레포 1(exercise)에서
# 검증된 WSL2 venv 를 재사용하고, 돌리기 전에 환경을 검사합니다.
# 스캔 자체는 수 분이지만 환경이 틀리면 import 에서 죽거나 — 더 나쁘게 —
# 사슬을 잘못 잡아 조용히 틀린 ddG 를 냅니다. 그래서 먼저 검사합니다.
#
# 사용
# ----
#   wsl                                   # Windows 터미널에서 WSL 진입
#   cd /mnt/c/Users/PC/Downloads/pdl1-binder-design
#   bash scripts/run_ala_scan.sh          # 검사 + 본실행
#   bash scripts/run_ala_scan.sh --check  # 검사만
#
# --check 뒤의 인자는 ala_scan.py 로 넘어갑니다. 비용을 먼저 재십시오:
#   bash scripts/run_ala_scan.sh --relax sc --residues A123,A32 --repeats 1
#
# venv 위치를 바꾸려면:  VENV=/경로/.venv-linux bash scripts/run_ala_scan.sh
set -euo pipefail

VENV="${VENV:-/mnt/c/Users/PC/Downloads/exercise/.venv-linux}"
PDB="${PDB:-inputs/PDL1_PD1_4ZQK.pdb}"
HOTSPOTS="${HOTSPOTS:-inputs/hotspots.json}"
REPEATS="${REPEATS:-3}"
OUT="${OUT:-ala_scan.csv}"
CHECK_ONLY=0
if [ "${1:-}" = "--check" ]; then CHECK_ONLY=1; shift; fi
# 남은 인자는 ala_scan.py 로 그대로 넘깁니다 (--relax sc 등)

die() { printf '\n[중단] %s\n' "$1" >&2; exit 1; }
ok()  { printf '  [ok] %s\n' "$1"; }

echo "=== 환경 검사 ==="

# 1) WSL 인가. Windows 쪽 python 으로 돌면 PyRosetta import 가 실패합니다.
grep -qi microsoft /proc/version 2>/dev/null || die "WSL2 안에서 실행하십시오. Windows 터미널에서 먼저 'wsl' 을 치십시오."
ok "WSL2"

# 2) venv 가 있는가. 없으면 만들지 않고 멈춥니다 — PyRosetta 설치는
#    라이선스 동의가 필요하고, 스크립트가 대신 할 일이 아닙니다.
[ -x "$VENV/bin/python" ] || die "venv 가 없습니다: $VENV
  레포 1의 venv 를 쓰거나, 아래 '처음 설치' 절차를 따르십시오."
PY="$VENV/bin/python"
ok "venv $("$PY" -V 2>&1)"

# 3) conda base 가 섞여 있지 않은가. 레포 1에서 겪은 함정입니다.
if [ -n "${CONDA_DEFAULT_ENV:-}" ]; then
  printf '  [주의] conda 환경 %s 이 활성입니다. venv 의 python 을 직접 쓰므로\n' "$CONDA_DEFAULT_ENV"
  printf '         보통 무해하지만, import 오류가 나면 conda deactivate 후 다시 하십시오.\n'
fi

# 4) PyRosetta 가 실제로 import 되는가. 버전까지 찍습니다.
"$PY" -c 'import pyrosetta,sys; print("  [ok] PyRosetta", pyrosetta.version() if hasattr(pyrosetta,"version") else "(버전 미확인)")' \
  2>/dev/null || die "PyRosetta import 실패. venv 가 맞는지 확인하십시오: $VENV"

# 5) 입력이 있는가. 그리고 사슬이 2개인가 — 여기서 걸러야
#    스크립트가 20분 돌고 나서 죽는 일이 없습니다.
for f in "$PDB" "$HOTSPOTS"; do
  [ -f "$f" ] || die "입력이 없습니다: $f"
done
ok "입력 $PDB / $HOTSPOTS"
# 명령 치환으로 받습니다. 파이프 + read 로 받으면 tr 의 출력에 개선문자가
# 없어 read 가 변수를 설정하고도 EOF 로 비-0 을 반환하고, set -e 가
# 거기서 스크립트를 죽입니다 (실제로 겪음).
chains=$(awk '$1=="ATOM"{print substr($0,22,1)}' "$PDB" | sort -u | tr -d '\n')
echo "  [ok] 사슬: $chains"
if [ "$chains" != "AB" ]; then
  printf '  [주의] 사슬이 A,B 가 아닙니다. ala_scan.py 는 A_B 계면을 가정합니다.\n'
fi
nA=$(awk '$1=="ATOM" && substr($0,22,1)=="A"' "$PDB" | wc -l)
nB=$(awk '$1=="ATOM" && substr($0,22,1)=="B"' "$PDB" | wc -l)
echo "  [ok] 원자 수: A $nA / B $nB"

echo
if [ "$CHECK_ONLY" = 1 ]; then echo "검사만 수행했습니다. 본실행: bash scripts/run_ala_scan.sh"; exit 0; fi

echo "=== 본실행 ==="
echo "  반복 $REPEATS, 출력 $OUT"
echo "  22잔기 x ($REPEATS+1) 회 계면 재계산입니다. 수 분 걸립니다."
echo
# -I: 격리 모드. 작업 디렉터리의 동명 모듈을 import 하지 않습니다.
exec "$PY" -I scripts/ala_scan.py \
  --pdb "$PDB" --hotspots "$HOTSPOTS" --repeats "$REPEATS" --out "$OUT" "$@"
