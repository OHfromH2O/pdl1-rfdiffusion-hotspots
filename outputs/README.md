# outputs/ — 무엇이 커밋되어 있고 무엇이 아닌가

## 커밋된 것

| | 수 | 내용 |
|---|---|---|
| `pdl1_{epi,epi2,ctrl,ctrl2}_*.pdb` | 100 | **RFdiffusion 백본.** 백본만 있습니다(N/CA/C/O). 이것이 생성 단계의 1차 산출물입니다 |
| `pdl1_*.trb` | 100 | 각 백본의 번호 매핑. H2a 의 타깃 잔기 번호 복원에 **필수**입니다 |
| `h1_af2/epi_af2.sc`, `h1_af2/ctrl_af2.sc` | 2 | AF2 initial guess 의 원시 점수 — `pae_interaction`, `plddt_*`, RMSD |
| `h1_af2/{epi,ctrl}/pdl1_{epi_23,ctrl_12}_*.pdb` | 2 | 그림에 쓴 복합체. 클론만으로 `scripts/07_render.py` 가 돌도록 |

## 커밋되지 않은 것

**AF2 initial guess 가 낸 복합체 98건**(약 26 MB)입니다.

빼는 이유는 용량 하나뿐이고, **측정값을 숨기는 것이 아닙니다.** 이 구조들에서
잰 값은 전부 레포에 있습니다.

| 측정 | 어디에 |
|---|---|
| `pae_interaction`, `plddt_binder`, RMSD | `outputs/h1_af2/*_af2.sc` |
| `dG_separated`, `hbonds_int`, `dSASA_int` | `../h1_rosetta.csv` |
| 위 둘을 합친 표 | `../h1_joined.csv` |
| 핫스팟 접촉 (H2b, 중원자 4.5 Å) | `../h2_af2_heavy.csv` |
| 핫스팟 접촉 (CA–CA 8 Å, 보조) | `../h2_af2_ca8.csv` |

## 재생성

커밋된 RFdiffusion 백본 100건에서 시작해
[`../notebooks/H1_dl_binder_design.ipynb`](../notebooks/H1_dl_binder_design.ipynb)
를 돌리면 복합체 100건이 다시 만들어집니다. 실측 비용은 Colab L4 기준
**ProteinMPNN 병렬 25.7분 + AF2 군별 59분 × 2**입니다. 설치에서 겪은 15개
실패 지점은 [`../notes/dl_binder_design_setup.md`](../notes/dl_binder_design_setup.md)
에 정리되어 있습니다.

생성 단계(RFdiffusion)부터 다시 하려면 Colab T4 에서 **1.8분/설계 × 100 ≈
3시간**이고, 설정은 [`../preregistration.md`](../preregistration.md) §2 에
고정되어 있습니다. 다만 확산 모델이라 **같은 백본이 다시 나오지는 않습니다.**
결론을 재현하려면 커밋된 백본에서 시작하십시오.
