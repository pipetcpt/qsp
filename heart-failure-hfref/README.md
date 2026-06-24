# 심부전 (감소 박출률, HFrEF) (Heart Failure with Reduced EF, HFrEF) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 심혈관

[![HFrEF QSP Model](hfref_qsp_model.png)](hfref_qsp_model.svg)

## 개요 (Overview)

감소 박출률 심부전(HFrEF, EF < 40%)은 신경호르몬(RAAS·교감신경계) 과활성이 심실 리모델링을 가속화하는 악순환 구조가 핵심 병태생리입니다. 전 세계 약 2,600만 명이 이환되어 있으며, 5년 사망률은 약 50%에 달합니다. 심근 손상 후 활성화된 AngII·알도스테론·노르에피네프린은 초기에 심박출량을 유지하지만, 만성화 시 심실 확장·섬유화·비대를 유발합니다. 현재 ARNI(사쿠비트릴-발사르탄)·베타차단제·MRA·SGLT2i의 4대 근거 기반 치료(GDMT)가 생존율을 현저히 개선합니다. 심박수 조절(이바브라딘)과 심실보조장치(VAD)가 추가 치료 선택지입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| RAAS 활성화 | 신장 저관류 → 레닌 → AngI → ACE → AngII → AT1R → 알도스테론 | 나트륨 저류·혈관 수축·심근 비대 |
| 교감신경계 항진 | 압수용체 반응↓ → NE 분비↑ → β₁R → 심박수·수축력↑ | 만성 심실 과부하, 베타1수용체 하향조절 |
| 나트륨이뇨펩티드 보상 | LVEDP↑ → BNP·ANP 분비 → NEP 분해 → 이뇨·혈관확장 | NT-proBNP 상승 |
| cGMP-PKG 신호 | ANP/BNP → pGC, NO → sGC → cGMP → PKG → 심근 이완 | ARNI로 cGMP 증폭 |
| 심실 리모델링 | TGF-β→ 섬유아세포→ ECM→ 심실 섬유화·확장 | LVEF 저하, 심실 확대 |
| 염증 활성화 | TNF-α·IL-6 → 심근 세포사멸·심실 기능 악화 | 심기능 저하 가속 |

## 주요 약물 표적 (Drug Targets)

- **ARNI(사쿠비트릴-발사르탄)**: NEP 억제(BNP 분해↓ → cGMP↑) + AT1R 차단; PARADIGM-HF에서 에나프릴 대비 사망·입원 감소
- **베타차단제(카르베딜롤·메토프롤롤)**: β₁/β₂R 차단 → 심박수·에너지 소비 감소 → 역리모델링; COPERNICUS·MERIT-HF 근거
- **MRA(에플레레논·스피로노락톤)**: 알도스테론 수용체 차단 → 항섬유화·이뇨; RALES·EMPHASIS-HF 근거
- **SGLT2 억제제(다파글리플로진·엠파글리플로진)**: 삼투 이뇨·심근 에너지 최적화·심실 부하 감소; DAPA-HF·EMPEROR-Reduced 근거
- **이바브라딘(Ivabradine)**: HCN(If 채널) 차단 → 심박수 선택적 감소; SHIFT 시험 근거
- **이뇨제(푸로세마이드)**: 과부하 증상 완화

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [hfref_qsp_model.dot](hfref_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 9 클러스터) |
| [hfref_qsp_model.svg](hfref_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [hfref_qsp_model.png](hfref_qsp_model.png) | PNG 이미지 (150 dpi) |
| [hfref_mrgsolve_model.R](hfref_mrgsolve_model.R) | mrgsolve ODE 모델 (약 26 구획 / 5개 치료 시나리오 + 용량-반응 분석) |
| [hfref_shiny_app.R](hfref_shiny_app.R) | Shiny 대시보드 |
| [hfref_references.md](hfref_references.md) | 참고문헌 (약 62편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: RAAS(AngI·AngII·Ang1-7·알도스테론), 교감신경(NE), BNP·NT-proBNP·cGMP, 혈역학(LVEDV·HR·SVR·LVEF), 리모델링(TGF-β1·섬유화·비대), 염증(TNF-α·IL-6), 약물 PK(LBQ657·발사르탄·베타차단제·MRA·SGLT2i·이바브라딘) 구획 포함(총 26 구획)
- **주요 치료 시나리오**: ① 무치료 기저선, ② ACEi+베타차단제(구 표준), ③ ARNI+BB+MRA(3제), ④ ARNI+BB+MRA+SGLT2i(4대 GDMT), ⑤ 최대 GDMT+이바브라딘(5제)
- **보정/근거**: PARADIGM-HF(McMurray 2014), DAPA-HF(McMurray 2019), EMPEROR-Reduced(Packer 2020), SHIFT(Swedberg 2010) 임상시험 데이터를 기반으로 LVEF 회복·NT-proBNP 감소·사망률을 정성적으로 보정

## Shiny 대시보드 (Dashboard)

환자 프로파일(기저 LVEF, BNP, NYHA 등급, 원인 선택) 탭, 약물 PK 동역학(ARNI·SGLT2i·BB 농도), RAAS·교감신경계 PD 지표, 심장 기능 임상 엔드포인트(LVEF·NT-proBNP·CO·HR), 5개 치료 시나리오 비교, 용량-반응 분석, 바이오마커(NT-proBNP·BNP·MAP) 탭으로 구성됩니다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("hfref_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("hfref_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg hfref_qsp_model.dot -o hfref_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [hfref_references.md](hfref_references.md) 참조 (약 62편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
