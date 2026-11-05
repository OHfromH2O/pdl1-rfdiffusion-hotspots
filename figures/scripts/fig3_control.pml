# 자동 생성 — scripts/07_render.py. 직접 수정하지 마십시오.
# .pml 에는 레포 기준 상대 경로만 씁니다. 절대 경로를 박으면 커밋된
# 파일이 다른 사람의 기계에서 못 돕니다. **레포 루트에서 실행하십시오**
# (scripts/07_render.py 로 돌리면 자동으로 그렇게 됩니다).
reinitialize
set ray_opaque_background, 0
set ray_shadows, 0
set antialias, 2
set cartoon_fancy_helices, 1
set field_of_view, 20.0
# auto_zoom 을 끕니다 — 켜 두면 load 마다 카메라가 리셋되어 뒤에서 지정한
# set_view 가 중간 load 에 덮입니다 (레포 1 에서 겪은 함정).
set auto_zoom, off
set ray_interior_color, grey80
set two_sided_lighting, on
bg_color white

# Fig 3 — 대조군 1위 pdl1_ctrl_12 (비-에피토프 핫스팟 지정)
#   dG_separated -20.25 REU   dSASA_int 1281 A^2
#   pae_interaction 19.77   핫스팟 접촉 0/4
#
# 번호 보정: 이 PDB 의 사슬 B 는 86- 에서 시작합니다(binder 뒤로
# 밀림). 오프셋 68 를 빼서 PD-L1 원본 번호로 되돌립니다.
load inputs/PDL1_PD1_4ZQK.pdb, nat
load outputs/h1_af2/ctrl/pdl1_ctrl_12_dldesign_0_cycle1_af2pred.pdb, des

alter des and chain B, resi=str(int(resi)-68)
sort

# **native 를 설계에 맞춥니다.** 반대로 하면 설계 좌표가 변해 아래에서
# 계산한 set_view 가 어긋납니다.
align nat and chain A and name CA, des and chain B and name CA

hide everything

show surface, des and chain B
color grey80, des and chain B
color red, des and chain B and resi 32+106+109+132
color marine, des and chain B and resi 56+115+121+123

# 설계한 binder
show cartoon, des and chain A
color orange, des and chain A

# native PD-1 — 실제 에피토프의 위치 표시
show cartoon, nat and chain B
color palecyan, nat and chain B
set cartoon_transparency, 0.7, nat and chain B

# 지정한 핫스팟 면이 정면에 오도록 계산한 카메라. 복합체 전체를 담습니다.
# 카메라 검증 (화면 좌표, 0 = 중앙). 패치 중심은 반드시 (0.00, 0.00):
#   패치 중심: 가로 -0.00  세로 -0.00
#   binder 중심: 가로 -0.06  세로 +0.07
set_view (\
     -0.2264,     0.9740,     0.0000,\
      0.2497,     0.0580,     0.9666,\
      0.9415,     0.2188,    -0.2564,\
      0.0000,     0.0000,  -311.5346,\
      2.6442,     1.7205,     0.9589,\
    279.6131,   349.2071,     0.0000 )

png figures/fig3_control.png, width=1800, height=1350, dpi=300, ray=1
