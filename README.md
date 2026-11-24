# C3 및 C4 식물의 rbcL 비교 분석

30종의 식물에서 엽록체 rbcL의 염기·코돈 조성과 Rubisco 큰 소단위체의 아미노산 서열을 비교하는 대학 생명과학 연구 과제입니다. C3·C4 각각 15종을 외떡잎식물 8종과 진정쌍떡잎식물 7종으로 구성합니다.

| 주차 | 연구 단계 |
| --- | --- |
| [11주차](week11/README.md) | 광합성 학습, 문헌 조사, 후보종 선정과 연구 설계 |
| [12주차](week12/README.md) | 네 그룹에서 1종씩 시험 수집 후 전체 데이터 검증과 GC 분석 |
| [13주차](week13/README.md) | 아미노산 조성, 단백질 특성, 코돈 사용 및 RSCU |
| [14주차](week14/README.md) | 다중서열정렬, 잔기 비교 및 rbcL 계통수 |
| [15주차](week15/README.md) | 그룹별 통계, 계통학적 한계 검토 및 최종 보고서 |

[연구 계획](week11/research_plan.md) · [후보종 CSV](data/processed/species_metadata.csv) · [진행 기록](PROGRESS.md)

## 실행 환경

Windows Anaconda PowerShell에서 `conda env create -f environment.yml`로 환경을 재현한 후 `conda activate rbcL`과 `jupyter lab`을 실행합니다. 기존에 rbcL 환경이 있으면 활성화해 사용합니다. Python 패키지 버전은 `requirements.txt`와 `environment.yml`에 기록합니다. 후보 메타데이터 검증은 프로젝트 루트에서 `python -m unittest discover -s tests -v`로 실행합니다.

14주차 정렬에는 Windows용 MUSCLE 5.3을 사용합니다. 공식 배포 파일을 프로젝트 루트의 `.tools/muscle.exe`에 두고 Windows PowerShell에서 `& .\.tools\muscle.exe -version`으로 확인합니다. 설치 파일의 SHA-256은 `a86637db99deb3efa1cb6762c4e3160ca33644608270eb167e4214f1f9409711`입니다. [공식 MUSCLE 5.3 배포](https://github.com/rcedgar/muscle/releases/tag/v5.3)

GenBank 원본은 `data/raw/`, 검증·처리 데이터는 `data/processed/`, 주차별 표와 그림은 각 주차의 `results/`와 `figures/`에 둡니다. 후보종의 Taxonomy ID, accession, 길이, 수집일과 NCBI 레코드 URL은 후속 수집에서 검증한 뒤 채웁니다.
