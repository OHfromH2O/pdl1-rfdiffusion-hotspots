#!/usr/bin/env python3
"""단계 7 — 구조 그림을 PyMOL 스크립트로 생성하고 헤드리스 렌더링한다.

레포 1의 `scripts/06_render.py`와 같은 방식입니다 — **수작업 스크린샷이
아니라 데이터에서 `.pml`을 만들어** 렌더링하고, 그 `.pml`을
`figures/scripts/`에 남겨 커밋합니다. 그림이 어떤 선택을 거쳐 만들어졌는지
제3자가 읽어서 검증할 수 있어야 하기 때문입니다.

대표 설계를 **손으로 고르지 않습니다.** `h1_joined.csv`에서 군별
`dG_separated` 최저값을 뽑습니다. 대조군 그림도 같은 규칙으로 고른
**대조군의 최고 설계**를 씁니다 — 그림이 결론을 꾸미지 않도록.

카메라를 계산합니다 (중요)
--------------------------
처음에는 `orient <핫스팟>`으로 카메라를 잡았는데 **빨강 패치가 보이지
않았습니다.** `orient`는 선택 영역의 주축을 화면축에 맞출 뿐, 그 면이
카메라를 **향하는지 등지는지**는 정하지 않습니다 — 반반입니다. 실제로
fig1·fig2 모두 패치가 뒷면으로 갔습니다.

그래서 추측을 버리고 **시점 행렬을 좌표에서 직접 계산해** `set_view`로
박습니다. 단백질 중심 → 패치 중심 방향을 카메라 축으로 삼으므로 패치가
**반드시** 정면에 옵니다. 결정론적이고, 다시 돌려도 같은 그림이 나옵니다.

중첩 방향도 뒤집었습니다 — 설계를 native에 맞추면 설계 좌표가 변해
계산한 카메라가 어긋납니다. **native를 설계에 맞춥니다.**

번호 보정
---------
AF2 출력은 `-maintain_res_numbering`을 주지 않아 타깃 사슬이 binder 뒤로
밀려 있습니다(실측 오프셋 53–83, 설계마다 다름). 각 PDB에서 오프셋을 읽어
`alter`로 **PD-L1 원본 번호(18–132)로 되돌립니다.** 그래야 선택문의
56·115·121·123이 보고서의 번호와 같은 것을 가리킵니다. 같은 번호를
잘못 읽는 결함이 `h2_hotspot_contacts.py`에 실제로 있었습니다(한계 20).

실행 (WSL2 — PyMOL은 네이티브 Windows에서 돌지 않습니다)
    python scripts/07_render.py
    python scripts/07_render.py --only 2
    python scripts/07_render.py --dry-run
"""
from __future__ import annotations

import argparse
import csv
import math
import pathlib
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
FIGDIR = ROOT / "figures"
PMLDIR = FIGDIR / "scripts"

NATIVE = ROOT / "inputs/PDL1_PD1_4ZQK.pdb"
NATIVE_REL = "inputs/PDL1_PD1_4ZQK.pdb"
JOINED = ROOT / "h1_joined.csv"
H2B = ROOT / "h2_af2_heavy.csv"
AF2DIR = ROOT / "outputs/h1_af2"

# 사전 등록된 두 핫스팟 묶음 (PD-L1 원본 번호)
EPI = [56, 115, 121, 123]
CTRL = [32, 106, 109, 132]

# 4ZQK: A = PD-L1 (18-132), B = PD-1 (33-146). 실측 확인했습니다.
NAT_TARGET, NAT_BINDER = "A", "B"

NL = chr(10)          # f-string 안에서 역슬래시를 쓰지 않기 위해
FOV = 20.0          # PyMOL 기본 field_of_view
WIDTH, HEIGHT = 1800, 1350

PREAMBLE = f"""# 자동 생성 — scripts/07_render.py. 직접 수정하지 마십시오.
# .pml 에는 레포 기준 상대 경로만 씁니다. 절대 경로를 박으면 커밋된
# 파일이 다른 사람의 기계에서 못 돕니다. **레포 루트에서 실행하십시오**
# (scripts/07_render.py 로 돌리면 자동으로 그렇게 됩니다).
reinitialize
set ray_opaque_background, 0
set ray_shadows, 0
set antialias, 2
set cartoon_fancy_helices, 1
set field_of_view, {FOV}
# auto_zoom 을 끕니다 — 켜 두면 load 마다 카메라가 리셋되어 뒤에서 지정한
# set_view 가 중간 load 에 덮입니다 (레포 1 에서 겪은 함정).
set auto_zoom, off
set ray_interior_color, grey80
set two_sided_lighting, on
bg_color white
"""


# ---------------------------------------------------------------------
# 좌표와 카메라
# ---------------------------------------------------------------------
def read_atoms(pdb: pathlib.Path, chain: str, resi: set[int] | None = None):
    """(x, y, z) 목록. resi 가 주어지면 그 잔기만."""
    out = []
    with open(pdb) as fh:
        for l in fh:
            if not l.startswith("ATOM") or l[21] != chain:
                continue
            if resi is not None and int(l[22:26]) not in resi:
                continue
            out.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return out


def _centroid(p):
    n = len(p)
    return (sum(a[0] for a in p) / n, sum(a[1] for a in p) / n,
            sum(a[2] for a in p) / n)


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _unit(a):
    n = math.sqrt(sum(x * x for x in a))
    if n < 1e-9:
        sys.exit("카메라 축이 0 입니다 — 패치 중심과 단백질 중심이 같습니다")
    return (a[0] / n, a[1] / n, a[2] / n)


def view_block(all_pts, patch_pts, center_pts=None, pad=1.06) -> str:
    """패치가 **반드시 정면에 오는** set_view 블록을 만든다.

    카메라의 z축(장면에서 카메라로 향하는 방향)을 '단백질 중심 → 패치
    중심'으로 잡습니다. 그러면 패치가 카메라 쪽 반구에 놓입니다.
    `orient`는 이 부호를 정해 주지 않아 패치가 뒷면으로 가곤 합니다.

    set_view 의 18개 값: 회전행렬 9(행 우선) + 카메라 좌표계 위치 3
    + 회전 원점(월드) 3 + near + far + orthoscopic.
    """
    # **두 역할을 분리합니다.**
    #   center_pts : 카메라 방향을 정하는 기준 (타깃 단백질만)
    #   all_pts    : 화면에 담겨야 하는 것 전부 (참조 분자 포함)
    # 하나로 합치면 참조 분자가 중심을 끌어당겨 패치가 가장자리로
    # 밀립니다 — fig1b 에서 실제로 그랬습니다.
    c = _centroid(center_pts if center_pts else all_pts)
    z = _unit(_sub(_centroid(patch_pts), c))
    up = (0.0, 0.0, 1.0) if abs(z[2]) < 0.9 else (0.0, 1.0, 0.0)
    x = _unit(_cross(up, z))
    y = _cross(z, x)

    # 거리는 **화면에 투영된 범위**로 정합니다. 3차원 반지름을 쓰면
    # ctrl_12 처럼 길쭉한 복합체가 과도하게 축소되어 여백만 커집니다
    # (실제로 그렇게 나왔습니다). 각 원자를 카메라 좌표 (u, v, w) 로
    # 바꾸고, 세로·가로 화각 안에 들어오는 최소 거리를 구합니다.
    tv = math.tan(math.radians(FOV / 2))        # 세로 half-FOV
    th = tv * (WIDTH / HEIGHT)                  # 가로 half-FOV
    uvw = [(sum(a * b for a, b in zip(_sub(p, c), x)),
            sum(a * b for a, b in zip(_sub(p, c), y)),
            sum(a * b for a, b in zip(_sub(p, c), z))) for p in all_pts]
    dist = pad * max(max(w + abs(v) / tv, w + abs(u) / th) for u, v, w in uvw)
    depth = [dist - w for _, _, w in uvw]
    near = max(1.0, min(depth) - 3.0)
    far = max(depth) + 3.0

    m = (list(x) + list(y) + list(z) + [0.0, 0.0, -dist] + list(c)
         + [near, far, 0.0])
    rows = [", ".join(f"{v:10.4f}" for v in m[i:i + 3]) for i in range(0, 18, 3)]
    return "set_view (\\\n  " + ",\\\n  ".join(rows) + " )"


def check_view(blk: str, pts_named) -> list[str]:
    """생성한 set_view 로 점들을 화면 좌표(-1..1)에 투영해 본다.

    **눈으로 판정하지 않기 위해서입니다.** 작업 중 저는 그림을 보고
    "패치가 흩어져 있다", "패치가 치우쳐 있다"고 두 번 적었고 두 번 다
    재 보니 틀렸습니다. 렌더 결과의 옳고 그름은 이 수치로 판정합니다.

    패치 중심은 시야축 위에 있도록 만들었으므로 (0.00, 0.00) 이어야
    합니다. 그렇지 않으면 카메라 계산이 깨진 것입니다.
    """
    import re
    v = [float(x) for x in
         re.search(r"set_view \((.*?)\)", blk, re.S).group(1)
         .replace("\\", " ").split(",")]
    R, dist, c = [v[0:3], v[3:6], v[6:9]], -v[11], v[12:15]
    tv = math.tan(math.radians(FOV / 2))
    th = tv * (WIDTH / HEIGHT)
    out = []
    for name, p in pts_named:
        d = [p[i] - c[i] for i in range(3)]
        uu = sum(R[0][i] * d[i] for i in range(3))
        vv = sum(R[1][i] * d[i] for i in range(3))
        w = sum(R[2][i] * d[i] for i in range(3))
        z = dist - w
        out.append(f"{name}: 가로 {uu / (th * z):+.2f}  세로 {vv / (tv * z):+.2f}"
                   + ("  [화면 밖]" if max(abs(uu / (th * z)),
                                          abs(vv / (tv * z))) > 1.0 else ""))
    return out


def sel(nums: list[int]) -> str:
    return "+".join(str(n) for n in nums)


def target_offset(pdb: pathlib.Path) -> int:
    """타깃 사슬(B) 첫 잔기 번호 - 18. 설계마다 다릅니다."""
    nums = [int(l[22:26]) for l in open(pdb)
            if l.startswith("ATOM") and l[21] == "B"]
    if not nums:
        sys.exit(f"{pdb.name}: 사슬 B 가 없습니다")
    return min(nums) - 18


def pick(group: str) -> dict:
    """군별 dG_separated 최저 설계. 손으로 고르지 않습니다."""
    if not JOINED.is_file():
        sys.exit(f"{JOINED.name} 가 없습니다. 먼저 scripts/h1_join.py 를 돌리십시오.")
    rows = [r for r in csv.DictReader(open(JOINED, encoding="utf-8"))
            if r["group"] == group]
    if not rows:
        sys.exit(f"{group} 행이 없습니다")
    best = min(rows, key=lambda r: float(r["dG"]))
    hit = ""
    if H2B.is_file():
        m = {r["name"]: r for r in
             csv.DictReader(open(H2B, encoding="utf-8"))}.get(best["name"])
        if m:
            hit = f"{m['n_hotspot_hit']}/4"
    sub = "epi" if group == "epitope" else "ctrl"
    pdb = AF2DIR / sub / f"{best['name']}.pdb"
    if not pdb.is_file():
        sys.exit(f"{pdb} 가 없습니다")
    return dict(name=best["name"].replace("_dldesign_0_cycle1_af2pred", ""),
                pdb=pdb, rel=pdb.relative_to(ROOT).as_posix(),
                dG=float(best["dG"]), pae=float(best["pae"]),
                sasa=float(best["sasa"]), hit=hit, offset=target_offset(pdb))


# ---------------------------------------------------------------------
# Fig 1 — 두 패치가 서로 다른 면에 있다
# ---------------------------------------------------------------------
def pml_fig1():
    """두 장을 **각각 다른 .pml 로** 냅니다.

    처음에는 한 .pml 에서 set_view -> png 를 두 번 했는데, **두 번째
    그림만 잘렸습니다.** 수학상 두 시점 모두 세로 채움 0.94 로 동일한데
    렌더 결과가 달랐으므로 PyMOL 쪽 상태가 남는 것입니다. 원인을 더
    파고드는 대신 변수를 없앱니다 — 그림 하나에 .pml 하나.

    축척 계산에 **PD-1(사슬 B)도 넣습니다.** 사슬 A 만으로 맞췄더니
    fig1a 에서 PD-1 이 왼쪽으로 잘렸습니다.
    """
    tgt_only = read_atoms(NATIVE, NAT_TARGET)
    allp = tgt_only + read_atoms(NATIVE, NAT_BINDER)
    epi = read_atoms(NATIVE, NAT_TARGET, set(EPI))
    ctrl = read_atoms(NATIVE, NAT_TARGET, set(CTRL))
    if not epi or not ctrl:
        sys.exit("4ZQK 에서 핫스팟 잔기를 찾지 못했습니다")
    gap = math.dist(_centroid(epi), _centroid(ctrl))

    common = PREAMBLE + f"""
# Fig 1 — 4ZQK. PD-L1(사슬 A) 표면 위의 두 핫스팟 묶음.
#
# 답하는 질문: 대조군 패치가 **정말 다른 자리**인가.
# 두 패치 중심 사이 거리 {gap:.1f} A. 반대 면이라 한 장에 담기지 않아
# 앞/뒤 두 장으로 냅니다. 카메라는 계산된 set_view 이고 orient 가 아닙니다.
load {NATIVE_REL}, nat

hide everything
show surface, nat and chain {NAT_TARGET}
color grey80, nat and chain {NAT_TARGET}
color red, nat and chain {NAT_TARGET} and resi {sel(EPI)}
color marine, nat and chain {NAT_TARGET} and resi {sel(CTRL)}

# native PD-1 — 실제 에피토프의 위치 표시. 얇게 둡니다.
show cartoon, nat and chain {NAT_BINDER}
color palecyan, nat and chain {NAT_BINDER}
set cartoon_transparency, 0.55, nat and chain {NAT_BINDER}
"""
    out = []
    for key, patch, note in (
            ("fig1a_epitope_face", epi,
             "에피토프 면 — PD-1 이 닿는 자리가 정면에 옵니다"),
            ("fig1b_control_face", ctrl,
             "대조 면 — 같은 분자를 반대쪽에서. PD-1 은 뒤로 갑니다")):
        blk = view_block(allp, patch, tgt_only)
        chk = check_view(blk, [("패치 중심", _centroid(patch))]
                         + [(f"A{n}", _centroid(read_atoms(NATIVE, NAT_TARGET,
                                                           {n}))) for n in
                            (EPI if key.startswith("fig1a") else CTRL)])
        chk_txt = NL.join("#   " + t for t in chk)
        body = common + f"""
# {note}
# 카메라 검증 (화면 좌표, 0 = 중앙 / ±1 = 가장자리). 눈으로 보지 않고
# 이 수치로 판정합니다 — 패치 중심은 반드시 (0.00, 0.00) 이어야 합니다.
{chk_txt}
{blk}
png figures/{key}.png, width={WIDTH}, height={HEIGHT}, dpi=300, ray=1
"""
        out.append((key, body, [f"{key}.png"]))
    return out


# ---------------------------------------------------------------------
# Fig 2·3 — 대표 복합체. 두 그림의 구성을 같게 유지합니다.
# ---------------------------------------------------------------------
def pml_design(d: dict, hotspots: list[int], other: list[int],
               fname: str, title: str):
    off = d["offset"]
    # 카메라는 **설계 좌표**에서 계산합니다. 그래서 설계를 움직이지 않고
    # native 를 설계에 맞춥니다 (반대로 하면 카메라가 어긋납니다).
    tgt = read_atoms(d["pdb"], "B")
    bnd = read_atoms(d["pdb"], "A")
    patch = read_atoms(d["pdb"], "B", {n + off for n in hotspots})
    if not patch:
        sys.exit(f"{d['name']}: 핫스팟 잔기를 찾지 못했습니다 (오프셋 {off})")

    blk = view_block(tgt + bnd, patch, tgt)
    chk_txt = NL.join("#   " + t for t in check_view(
        blk, [("패치 중심", _centroid(patch)),
              ("binder 중심", _centroid(bnd))]))
    body = PREAMBLE + f"""
# {title}
#   dG_separated {d['dG']:+.2f} REU   dSASA_int {d['sasa']:.0f} A^2
#   pae_interaction {d['pae']:.2f}   핫스팟 접촉 {d['hit'] or 'n/a'}
#
# 번호 보정: 이 PDB 의 사슬 B 는 {off + 18}- 에서 시작합니다(binder 뒤로
# 밀림). 오프셋 {off} 를 빼서 PD-L1 원본 번호로 되돌립니다.
load {NATIVE_REL}, nat
load {d['rel']}, des

alter des and chain B, resi=str(int(resi)-{off})
sort

# **native 를 설계에 맞춥니다.** 반대로 하면 설계 좌표가 변해 아래에서
# 계산한 set_view 가 어긋납니다.
align nat and chain {NAT_TARGET} and name CA, des and chain B and name CA

hide everything

show surface, des and chain B
color grey80, des and chain B
color red, des and chain B and resi {sel(hotspots)}
color marine, des and chain B and resi {sel(other)}

# 설계한 binder
show cartoon, des and chain A
color orange, des and chain A

# native PD-1 — 실제 에피토프의 위치 표시
show cartoon, nat and chain {NAT_BINDER}
color palecyan, nat and chain {NAT_BINDER}
set cartoon_transparency, 0.7, nat and chain {NAT_BINDER}

# 지정한 핫스팟 면이 정면에 오도록 계산한 카메라. 복합체 전체를 담습니다.
# 카메라 검증 (화면 좌표, 0 = 중앙). 패치 중심은 반드시 (0.00, 0.00):
{chk_txt}
{blk}

png figures/{fname}, width={WIDTH}, height={HEIGHT}, dpi=300, ray=1
"""
    return body, [fname]


def run_pymol(pml: pathlib.Path) -> None:
    """헤드리스 실행. -c = GUI 없음, -q = 배너 억제."""
    exe = shutil.which("pymol")
    cmd = ([exe, "-cq", str(pml)] if exe
           else [sys.executable, "-m", "pymol", "-cq", str(pml)])
    # .pml 이 상대 경로를 쓰므로 레포 루트에서 띄웁니다.
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    if r.returncode != 0:
        print(r.stdout[-2000:], file=sys.stderr)
        print(r.stderr[-2000:], file=sys.stderr)
        raise SystemExit(f"PyMOL 실패: {pml.name}")
    for line in (r.stdout + r.stderr).splitlines():
        if "Error" in line or "error:" in line:
            print(f"  [warn] {line}")


def main() -> int:
    ap = argparse.ArgumentParser(description="단계 7 — 구조 그림")
    ap.add_argument("--only", type=int, choices=[1, 2, 3])
    ap.add_argument("--dry-run", action="store_true",
                    help=".pml 만 만들고 렌더링은 생략")
    a = ap.parse_args()

    FIGDIR.mkdir(exist_ok=True)
    PMLDIR.mkdir(parents=True, exist_ok=True)
    if not NATIVE.is_file():
        sys.exit(f"{NATIVE} 가 없습니다")

    jobs = []
    if a.only in (None, 1):
        jobs.extend(pml_fig1())
    if a.only in (None, 2):
        e = pick("epitope")
        print(f"[설계군 1위] {e['name']}  dG {e['dG']:+.2f} REU  "
              f"핫스팟 {e['hit'] or 'n/a'}  오프셋 {e['offset']}")
        body, pngs = pml_design(e, EPI, CTRL, "fig2_top_design.png",
                                f"Fig 2 — 설계군 1위 {e['name']} "
                                f"(에피토프 핫스팟 지정)")
        jobs.append(("fig2_top_design", body, pngs))
    if a.only in (None, 3):
        c = pick("control")
        print(f"[대조군 1위] {c['name']}  dG {c['dG']:+.2f} REU  "
              f"핫스팟 {c['hit'] or 'n/a'}  오프셋 {c['offset']}")
        body, pngs = pml_design(c, CTRL, EPI, "fig3_control.png",
                                f"Fig 3 — 대조군 1위 {c['name']} "
                                f"(비-에피토프 핫스팟 지정)")
        jobs.append(("fig3_control", body, pngs))

    for name, body, pngs in jobs:
        pml = PMLDIR / f"{name}.pml"
        pml.write_text(body, encoding="utf-8", newline="\n")
        print(f"  -> {pml.relative_to(ROOT)}")
        if a.dry_run:
            continue
        for p in pngs:        # 옛 그림을 지워 두어야 실패가 눈에 보입니다
            (FIGDIR / p).unlink(missing_ok=True)
        run_pymol(pml)
        for p in pngs:
            f = FIGDIR / p
            if f.is_file():
                print(f"  -> figures/{p} ({f.stat().st_size / 1e3:.0f} KB)")
            else:
                print(f"  [실패] {p} 이 만들어지지 않았습니다")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
