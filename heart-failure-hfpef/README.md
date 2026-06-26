# 심부전 (보존 박출률, HFpEF) (Heart Failure with Preserved EF, HFpEF) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 심혈관

[![HFpEF QSP Model](hfpef_qsp_model.png)](hfpef_qsp_model.svg)

## 개요 (Overview)

보존 박출률 심부전(HFpEF, EF ≥ 50%)은 전체 심부전 환자의 절반 이상을 차지하며, 고혈압·비만·당뇨·신장질환 등 대사 동반 질환에 의한 심근 경직(diastolic dysfunction), 전신 염증, 미세혈관 기능부전이 복합적으로 작용하는 확장기 부전 증후군입니다. Paulus-Tschöpe 패러다임에 따르면 대사 동반 질환 → 관상동맥 미세혈관 내피 염증 → cGMP-PKG 신호 저하 → 타이틴 경직 → 이완 장애의 경로가 핵심입니다. HFrEF와 달리 명확한 생존 이득이 입증된 치료제가 제한적이었으나, EMPEROR-Preserved 및 DELIVER 시험으로 SGLT2 억제제가 최초로 입원 감소 효과를 증명했습니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 대사 스트레스·전신 염증 | 비만·당뇨 → TNF-α·IL-6·CRP 상승 → 내피 산화 스트레스 | 미세혈관 기능부전 |
| cGMP-PKG 신호 저하 | eNOS↓·sGC↓ → cGMP↓ → PKG 활성↓ → 타이틴 인산화↓ | 심근 경직, LV 이완 장애 |
| RAAS 활성화 | Ang II → aldosterone → 심근 섬유화·나트륨 저류 | LV 비대, 부종 |
| 나트륨이뇨펩티드 과부하 | LVEDP 상승 → BNP·ANP 상승 → 보상 작동 | NT-proBNP 상승 |
| 타이틴 경직 | PKG 인산화 감소 → 티틴 N2B 과탄성 → 이완 불량 | E/e' 비 상승, 운동 불내성 |
| 심근 섬유화 | TGF-β→ 섬유아세포→ECM 침착 → 좌심실 경직 | 이완기 기능 악화 |

## 주요 약물 표적 (Drug Targets)

- **SGLT2 억제제(엠파글리플로진·다파글리플로진)**: 삼투 이뇨·혈역학 개선·심근 에너지 대사 최적화; EMPEROR-Preserved·DELIVER에서 심부전 입원 감소
- **MRA(스피로노락톤·에플레레논)**: 알도스테론 수용체 차단 → 심근 섬유화·나트륨 저류 억제; TOPCAT 일부 하위분석 근거
- **ARNI(사쿠비트릴-발사르탄)**: NEP 억제(BNP 분해 억제) + AT1R 차단 → cGMP↑; PARAGON-HF에서 경계 유의미 결과
- **이뇨제(푸로세마이드)**: 과부하 증상 완화; 직접 생존 이득 미입증
- **sGC 자극제(베리시구아트)**: sGC 직접 자극 → cGMP↑ → PKG↑ → 타이틴 연화; VITALITY-HFpEF 연구 중

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [hfpef_qsp_model.dot](hfpef_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 120+ 노드 / 14 클러스터) |
| [hfpef_qsp_model.svg](hfpef_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [hfpef_qsp_model.png](hfpef_qsp_model.png) | PNG 이미지 (150 dpi) |
| [hfpef_mrgsolve_model.R](hfpef_mrgsolve_model.R) | mrgsolve ODE 모델 (약 27 구획 / 5개 치료 시나리오) |
| [hfpef_shiny_app.R](hfpef_shiny_app.R) | Shiny 대시보드 |
| [hfpef_references.md](hfpef_references.md) | 참고문헌 (약 44편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK(엠파글리플로진 장·중심·말초, 사쿠비트릴-LBQ657 장·중심, 발사르탄 장·중심, 피네레논 장·중심, 푸로세마이드 장·중심) + 질환 PD(AngII·알도스테론 동태, ANP·BNP·cGMP·PKG, LV 비대지수·섬유화·타이틴 경직·LVEDP·SVR·IL-6·CRP·GFR·나트륨 배설·NT-proBNP) 구획 포함
- **주요 치료 시나리오**: ① 무치료 기저선, ② SGLT2i 단독(엠파글리플로진 10 mg), ③ SGLT2i+MRA 병합, ④ SGLT2i+ARNI 병합, ⑤ 4제 병합(SGLT2i+MRA+ARNI+이뇨제)
- **보정/근거**: EMPEROR-Preserved(Anker 2021, NEJM), DELIVER(Solomon 2022, NEJM), TOPCAT 하위분석, PARAGON-HF 데이터를 기반으로 NT-proBNP 감소 및 심부전 입원 감소율을 정성적으로 보정

## Shiny 대시보드 (Dashboard)

환자 프로파일(기저 LVEF, NT-proBNP, 동반 질환 선택) 탭, 약물 PK 동역학(SGLT2i·ARNI 농도), 심장 PD 지표(LVEDP·cGMP·타이틴 경직·섬유화), NT-proBNP·GFR·혈압 임상 엔드포인트, 심부전 입원 위험도, 5개 치료 시나리오 비교, 바이오마커(NT-proBNP·IL-6·CRP·eGFR) 탭으로 구성됩니다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("hfpef_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("hfpef_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg hfpef_qsp_model.dot -o hfpef_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [hfpef_references.md](hfpef_references.md) 참조 (약 44편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
